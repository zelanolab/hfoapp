% Modified from RippleLab
%   f_findHFOxEnEnt.m [As a part of HFO Detection Project]
%   Written by:
%   Miguel G. Navarrete Mejia
%   Electrical Engineering MS candidate
%   UNIVERSIDAD DE LOS ANDES
%   Colombia, 2012
%   mnavarretem@gmail.com


function evinfo = HFOAutoDetectorMNI( v, srate, cfg)
% 
%   v, row vector of time series   
%   srate, sampling rate, in Hz
%   cfg, structure contain following fields
%       * Fields accepted by HFO_FiltData
%       'bpfreq', band-pass frequency [lo, hi] Hz
%    	'bptype', filter type 'iir' | 'fir'
%       'notch', 'yes' | 'no', apply notch filter
%       'notch_freq', 50 | 60 line noise frequency
%       'notch_bandwidth', scalar, bandwidth of the filter
%       'notch_type, 'fir' | 'iir', notch filter type
%       * Fields for MNI detector
%        'epoch', scalar, epoch length in seconds
%        'hf_epcoh', scalar, Continous High Frequency Epoch, in seconds
%        'hf_percent', 0-1 scalar, Continous High Frequency Percentil Threshold     
%        'min_win', scalar, minimum HFO Time, in seconds   
%        'min_gap', scalar Minimum HFO Gap, in seconds
%        'thresh_percent', 0-1 scalar, Threshold precentil
%        'bs_win', scalar, Baseline window, in seconds
%        'bs_shift', scalar, Baseline Shift window in seconds
%        'bs_thresh', scalar, Baseline threshold
%        'bs_min', scalar, Baseline minimum time in seconds

evinfo = [];

%% Variable declarations
v_Freqs         = cfg.bpfreq;  % [lo, hi] band-pass frequency
s_Epoch         = cfg.epoch;                       % scalar, Cycle Time, in seconds
s_EpoCHF        = cfg.hf_epcoh;                          % Continous High Frequency Epoch, in seconds
s_PerCHF        = cfg.hf_percent;                      % Continous High Frequency Percentil Threshold       
s_MinWin        = cfg.min_win;                   % Minimum HFO Time  , in seconds   
s_MinGap        = cfg.min_gap;                   % Minimum HFO Gap, in seconds
s_ThresPerc     = cfg.thresh_percent;                   % Threshold precentil
s_BaseSeg       = cfg.bs_win;                    % Baseline window , in seconds
s_BaseShift     = cfg.bs_shift;                       % Baseline Shift window    
s_BaseThr       = cfg.bs_thresh;                         % Baseline threshold
s_BaseMin       = cfg.bs_min;	                        % Baseline minimum time, seconds
    
%% band-pass filter
v_SigFilt = HFOFiltData( v, srate, cfg);
            
%% RMS Calculus
s_Window        = round(0.002 * srate);
if mod(s_Window, 2) == 0
    s_Window = s_Window + 1;
end

v_Temp                      = v_SigFilt.^2;
v_Temp                      = filter(ones(1,s_Window),1, v_Temp)./s_Window;
v_RMS                       = zeros(numel(v_Temp), 1);
v_RMS(1:end - ceil(s_Window / 2) + 1) = v_Temp(ceil(s_Window / 2):end);
v_RMS                       = sqrt(v_RMS);
clear v_Temp

%% Baseline detection    
s_EpochSamples      = round(s_BaseSeg * srate);
s_FreqSeg           = numel(v_Freqs(1):5:v_Freqs(2));
s_StDevCycles       = 3;

s_WEMax = zeros(100,1); % Simulation of maximum baseline based in
                        % the calculus of uniform random values
for kk = 1 : 100 
    v_Segment       = rand(s_EpochSamples,1);
    v_AutoCorr      = xcorr(v_Segment)./sum(v_Segment.^2);
                            
	[m_WCoef] = HFOSpectrogram( v_AutoCorr, srate,...
        'foi', linspace( v_Freqs(1), v_Freqs(2), s_FreqSeg),...
        'ps_StDevCycles', s_StDevCycles);
    
    m_WCoef = flipud( abs(m_WCoef));                            
    m_ProbEnerg     = mean(m_WCoef.^2,2)/sum(mean(m_WCoef.^2,2));
    s_WEMax(kk)     = -sum(m_ProbEnerg.*log(m_ProbEnerg));
end

clear v_Segment m_ProbEnerg m_WCoef

s_WEMax             = median(s_WEMax);
v_InIndex           = 1:round(s_EpochSamples*s_BaseShift):numel(v);
v_EnIndex           = v_InIndex + s_EpochSamples - 1;
v_Idx               = v_EnIndex > numel(v);
v_InIndex(v_Idx)    = [];
v_EnIndex(v_Idx)    = [];
clear v_Idx

v_BaselineWindow	= zeros(size(v_SigFilt));

for kk = 1:numel(v_InIndex)              
    v_AutoCorr      = xcorr(v_SigFilt(v_InIndex(kk):v_EnIndex(kk)))./...
                    sum((v_SigFilt(v_InIndex(kk):v_EnIndex(kk))).^2);
                
	[m_WCoef] = HFOSpectrogram( v_AutoCorr, srate,...
        'foi', linspace( v_Freqs(1), v_Freqs(2), s_FreqSeg),...
        'ps_StDevCycles', s_StDevCycles);
    
    m_WCoef = flipud( abs( m_WCoef));
                            
    m_ProbEnerg 	= mean(m_WCoef.^2,2)/sum(mean(m_WCoef.^2,2));
    s_WESection     = -sum(m_ProbEnerg.*log(m_ProbEnerg));
    
    if s_WESection < s_BaseThr*s_WEMax
        v_BaselineWindow(v_InIndex(kk):v_EnIndex(kk))   = 1;
    end
end

clear m_ProbEnerg m_WCoef

%% Threshold calculation 
s_WindSamples	= (s_BaseMin/60) * numel(v);
s_MinWin        = s_MinWin * srate;
if sum(v_BaselineWindow) >= s_WindSamples
    v_Threshold     = f_PosBaseline();
else
    v_Threshold     = f_NegBaseline();
end

%% Event detection
v_EnergyThres   = v_RMS >= v_Threshold;

v_WindThres     = [0;v_EnergyThres;0];
v_WindJumps     = diff(v_WindThres);
v_WindJumUp     = find(v_WindJumps==1);
v_WindJumDown   = find(v_WindJumps==-1);
v_WindSelect    = find( (v_WindJumDown - v_WindJumUp) > s_MinWin);
if isempty(v_WindSelect)
    return
end

v_WindIni	= v_WindJumUp(v_WindSelect);
v_WindEnd	= v_WindJumDown(v_WindSelect);
s_MinGap	= s_MinGap * srate;

clear v_WindThres v_WindJumps v_WindJumUp v_WindJumDown
clear v_DistSelect v_WindSelect

while 1    
    if isempty(v_WindIni)
        return
    end
    
    if numel(v_WindIni) <2
        break
    end
    
    v_NextIni   = v_WindIni(2:end);
    v_LastEnd   = v_WindEnd(1:end-1);
    v_WinIdx	= (v_NextIni - v_LastEnd) < s_MinGap;
    
    if sum(v_WinIdx)==0
        break
    end
    v_NewEnd    = v_WindEnd(2:end);
    
    v_LastEnd(v_WinIdx) = v_NewEnd(v_WinIdx);
    v_WindEnd(1:end-1)  = v_LastEnd;
    
    v_Idx       = diff([0;v_WindEnd])~=0;
    v_WindIni   = v_WindIni(v_Idx);
    v_WindEnd   = v_WindEnd(v_Idx);
end

if ~isempty( v_WindIni)
    for k = 1 : length( v_WindIni)
        loc = v_WindIni( k) : v_WindEnd( k) - 1;
        evinfo.Location( k, :) = loc( [1, end]);
        [nb_cycle, avg_freq] = EstimateNbcycles( v_SigFilt( loc), srate);
        evinfo.NoCycles( k, 1) = nb_cycle;
        evinfo.PeakZScore( k, 1) = max( v_RMS( loc));
        evinfo.AvgFreq( k, 1) = avg_freq;
        evinfo.BPFreq( k, :) = cfg.bpfreq;
    end
end

clear v_WindSelect v_WindIni v_WindEnd v_WinDist

        
%% Functions      
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function v_Threshold = f_PosBaseline()     
        
        s_WindowThreshold   = round(s_Epoch * srate);
        v_BaselineWindow    = v_BaselineWindow(:);
        v_WindBaseline      = diff([0;v_BaselineWindow;0]);
        v_WindBaselineUp    = find(v_WindBaseline==1);         
        v_WindBaselineDown	= find(v_WindBaseline==-1)-1;
        
        v_IdxIni            = [];
        
        for jj = 1:numel(v_WindBaselineUp) 
            v_IdxAdd    = v_WindBaselineUp(jj):s_WindowThreshold:...
                        v_WindBaselineDown(jj);
                
            v_IdxIni    = horzcat(v_IdxIni,v_IdxAdd); %#ok<AGROW>
        end
        
        v_IdxEnd        = v_IdxIni + s_WindowThreshold -1;  
        
        v_IdxRem            = v_IdxEnd > numel(v_RMS);
        v_IdxIni(v_IdxRem)  = [];
        v_IdxEnd(v_IdxRem)  = [];

        v_Threshold     = zeros(size(v_RMS));
        
        for jj = 1:numel(v_IdxIni) 

            v_Section   = sort(v_RMS(v_IdxIni(jj):v_IdxEnd(jj)),'ascend');
            v_GamParams = fb_gamfit(v_Section);
            v_Percent   = fb_gamcdf(v_Section,v_GamParams(1),v_GamParams(2)); 
            s_Index     = find(v_Percent <= s_ThresPerc,1,'last');   
            v_Threshold(v_IdxIni(jj):v_IdxEnd(jj)) = v_Section(s_Index);
        end
        
    end
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function v_Threshold = f_NegBaseline()        
                            
        s_WindowCHF         = round(s_EpoCHF * srate);            
        
        v_IdxIni            = 1:s_WindowCHF:numel(v_RMS);
        v_IdxEnd            = v_IdxIni + s_WindowCHF - 1;
        v_IdxRem            = v_IdxEnd > numel(v_RMS);
        v_IdxIni(v_IdxRem)	= [];
        v_IdxEnd(v_IdxRem)	= [];
        
        v_Threshold     = zeros(size(v_RMS));
        
        for jj = 1:numel(v_IdxIni) 

            v_Section   = sort(v_RMS(v_IdxIni(jj):v_IdxEnd(jj)),'ascend');
            s_ThresLast = max(v_Section);
            
            while 1
                
                if sum(abs(v_Section)) == 0
                    break
                end
                
                v_GamParams = gamfit(v_Section);               
                v_Percent   = gamcdf(v_Section,v_GamParams(1),v_GamParams(2));
                s_ThresNew	= v_Section(find(v_Percent <= s_PerCHF,1,'last'));
                
                v_EnergyOver   = v_Section >= s_ThresNew;
                v_WindThr       = [0;v_EnergyOver;0];
                v_Jumps         = diff(v_WindThr);
                v_JumUp         = find(v_Jumps==1);
                v_JumDown       = find(v_Jumps==-1)-1;
                v_Dist          = v_JumDown - v_JumUp;
                
                v_Select        = (v_Dist > s_MinWin);                
                v_Select        = find(v_Select);
                
                if isempty(v_Select)
                    break
                end
                
                v_Ini       = v_JumUp(v_Select);
                v_End       = v_JumDown(v_Select);
                
                for ii = 1:numel(v_Ini)
                    v_Section(v_Ini(ii):v_End(ii))  = 0;
                    v_Section                       = sort(v_Section,'ascend');
                end
                
                s_ThresLast	= s_ThresNew;
                
            end  
            
            v_Threshold(v_IdxIni(jj):v_IdxEnd(jj)) = s_ThresLast;
        end
         
    end

end % main function