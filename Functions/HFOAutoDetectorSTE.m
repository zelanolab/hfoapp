% Modified from RippleLab
%   f_findHFOxSTF.m [As a part of HFO Detection Project]
%   Written by:
%   Miguel G. Navarrete Mejia
%   Electrical Engineering MS candidate
%   UNIVERSIDAD DE LOS ANDES
%   Colombia, 2012
%   mnavarretem@gmail.com


function evinfo = HFOAutoDetectorSTE( v, srate, cfg)
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
%       * Fields for STE detector
%        'rs_win', scalar, RMS window time in seconds
%        'rs_thresh', Threshold for RMS in standard deviation
%        'rs_min_dur',  Min window time for an HFO in seconds
%        'min_dist', Min Distance time Betwen two HFO candidates
%        'nocycle_thresh', % Minimum oscillations per interval
%        'peak_thresh', % Threshold for finding peaks
%        'epoch', % Cycle Time in seconds


v_Freqs         = cfg.bpfreq; % Filter freqs
s_Window        = cfg.rs_win ;           % RMS window time (ms)
s_RMSThresSD    = cfg.rs_thresh;                   % Threshold for RMS in standard deviation
s_MinWind       = cfg.rs_min_dur;             % Min window time for an HFO (ms)
s_MinTime       = cfg.min_dist;             % Min Distance time Betwen two HFO candidates
s_NumOscMin     = cfg.nocycle_thresh;                  % Minimum oscillations per interval
s_BPThresh      = cfg.peak_thresh;                   % Threshold for finding peaks
s_EpochLength   = cfg.epoch;                  % Cycle Time


%% Preprocessing Filter
v_SigFilt       = HFOFiltData(v, srate, cfg);

%% RMS Calculus
s_Window        = round(s_Window * srate);
if mod(s_Window, 2) == 0
    s_Window = s_Window + 1;
end
v_Temp                      = v_SigFilt.^2;
v_Temp                      = filter(ones(1,s_Window),1,v_Temp)./s_Window;
v_RMS                       = zeros(numel(v_Temp), 1);
v_RMS(1:end - ceil(s_Window / 2) + 1) = v_Temp(ceil(s_Window / 2):end);
v_RMS                       = sqrt(v_RMS);

clear v_Temp

%% Thresholding
s_MinWind       = round(s_MinWind * srate);
s_MinTime       = round(s_MinTime * srate);
s_EpochLength   = round(s_EpochLength * srate);
v_EpochTemp     = (1:s_EpochLength:length(v))';

if v_EpochTemp(end) < length(v)
    v_EpochTemp(end+1)  = length(v);
end

m_EpochLims     = [v_EpochTemp(1:end-1) v_EpochTemp(2:end)-1];
s_Epochs        = size(m_EpochLims,1);

clear v_EpochTemp s_EpochLength

m_HFOEvents = [];
s_BarCount  = 0;

for ii = 1:size(m_EpochLims,1)
    
    s_BarCount      = s_BarCount + 1;
    
    v_Window        = zeros(numel(v_RMS),1);
    v_Window(m_EpochLims(ii,1):m_EpochLims(ii,2)) = 1;
    
    v_RMSEpoch      = v_RMS.*v_Window;
    v_RMSInterval   = v_RMS(m_EpochLims(ii,1):m_EpochLims(ii,2));
    v_EpochFilt     = v_SigFilt(m_EpochLims(ii,1):m_EpochLims(ii,2));
    
    v_RMSThres      = v_RMSEpoch > ...
        (mean(v_RMSInterval)+ ...
        s_RMSThresSD*std(v_RMSInterval));
    
    if isempty(numel(find(v_RMSThres)))
        s_BarCount      = s_BarCount + 1;
        continue
    end
    
    v_WindThres     = [0;v_RMSThres;0];
    v_WindJumps     = diff(v_WindThres);
    v_WindJumUp     = find(v_WindJumps==1);
    v_WindJumDown   = find(v_WindJumps==-1);
    v_WinDist       = v_WindJumDown - v_WindJumUp;
    
    v_WindIni       = v_WindJumUp(v_WinDist > s_MinWind);
    v_WindEnd       = v_WindJumDown(v_WinDist > s_MinWind)-1;
    
    if isempty(v_WindIni)
        s_BarCount      = s_BarCount + 1;
        continue
    end
    
    clear v_WindThres v_WindJumps v_WindJumUp v_WindJumDown
    clear v_DistSelect v_WindSelect
    
    while 1
        v_NextIni   = v_WindIni(2:end);
        v_LastEnd   = v_WindEnd(1:end-1);
        v_WinIdx	= (v_NextIni - v_LastEnd) < s_MinTime;
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
    
    m_WindIntervals = [v_WindIni v_WindEnd];
    
    clear v_WindSelect v_WindIni v_WindEnd v_WinDist
    
    s_Count             = 1;
    m_WindSelect        = zeros(size(m_WindIntervals));
    
    s_Threshold         = mean(abs(v_EpochFilt)) + ...
        s_BPThresh.*std(abs(v_EpochFilt));
    s_TotalWindInterv	= size(m_WindIntervals,1);
    
    
    for jj=1:s_TotalWindInterv
        
        v_Temp          = abs(v_SigFilt(m_WindIntervals(jj,1):...
            m_WindIntervals(jj,2)));
        
        if numel(v_Temp) < 3
            continue
        end
        
        s_NumPeaks      = findpeaks(v_Temp,'minpeakheight',s_Threshold);
        clear v_Temp
        
        if isempty(s_NumPeaks) || length(s_NumPeaks) < s_NumOscMin
            continue;
        end
        
        m_WindSelect(s_Count,:) = [m_WindIntervals(jj,1)...
            m_WindIntervals(jj,2)];
        s_Count                 = s_Count + 1;
        
    end
    
    if any(m_WindSelect(:))
        m_HFOEvents     = vertcat(m_HFOEvents,...
            m_WindSelect(1:s_Count-1,:)); %#ok<AGROW>
    end
    
end


evinfo = [];
if ~isempty( m_HFOEvents)
    for k = 1 : size( m_HFOEvents, 1)
        loc = m_HFOEvents( k, 1) : m_HFOEvents( k, 2);
        evinfo.Location( k, :) = loc( [1, end]);
        [nb_cycle, avg_freq] = EstimateNbcycles( v_SigFilt( loc), srate);
        evinfo.NoCycles( k, 1) = nb_cycle;
        evinfo.PeakZScore( k, 1) = nan;
        evinfo.AvgFreq( k, 1) = avg_freq;
        evinfo.BPFreq( k, :) = cfg.bpfreq;
    end
end


end