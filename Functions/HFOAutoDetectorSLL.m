% Modified from RippleLab
%   f_findHFOxSLL.m [As a part of HFO Detection Project]
%   Written by:
%   Miguel G. Navarrete Mejia
%   Electrical Engineering MS candidate
%   UNIVERSIDAD DE LOS ANDES
%   Colombia, 2012
%   mnavarretem@gmail.com



function evinfo = HFOAutoDetectorSLL( v, srate, cfg)
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
%       * Fields for SLL detector
%        'win_len', scalar, length in seconds
%        'epoch', scalar, in seconds
%        'percent', 0-1, scalar
%        'min_dur', minimal duration, in seconds
%        'do_equal', 'yes' | 'no' Equalization Filter

v_Freqs = cfg.bpfreq;
s_Window        = round(cfg.win_len * srate);
s_EpochWind     = cfg.epoch;
s_Percentil     = cfg.percent;
s_MinWind       = cfg.min_dur; 


%% Equalization Filter and band-pass filter
if strcmpi( cfg.do_equal, 'yes')
    v_FiltSignal    = f_EqualizerFreqFilter(v,srate,1,v_Freqs);
    v_FiltSignal = HFOFiltData( v_FiltSignal, srate, cfg);

else
    v_FiltSignal = HFOFiltData( v, srate, cfg);
end

%% Energy Line Length
   
v_FiltSignal    = v_FiltSignal(:);
v_FiltSignal    = [v_FiltSignal(1);v_FiltSignal];
v_Temp          = abs(diff(v_FiltSignal));
v_Temp          = filter(ones(1,s_Window),1,v_Temp)./s_Window;
v_Energy        = zeros(numel(v_Temp), 1);
v_Energy(1:end - ceil(s_Window / 2) + 1) = 10*v_Temp(ceil(s_Window / 2):end);

% Window to the 1% for edge cutting
v_EdgeWind      = gausswin(round(numel(v_Energy)*.01),5.25);
v_EdgeWind      = v_EdgeWind(:);
[~,s_CenterIdx] = max(v_EdgeWind);
v_EdgeWind      = vertcat(v_EdgeWind(1:s_CenterIdx),...
                ones(numel(v_Energy)-numel(v_EdgeWind),1),...
                v_EdgeWind(s_CenterIdx+1:end));

v_Energy        = v_Energy.*v_EdgeWind;

           
%% Thresholding
  
s_MinWind       = round(s_MinWind * srate);
s_EpochWind     = round(s_EpochWind .* srate);
s_Epochs        = round(numel(v_Energy)/s_EpochWind);
v_Threshold     = zeros(numel(v_Energy),1);

for kk=1:s_Epochs

    s_Ini	= floor((kk-1) .* s_EpochWind)+1;
    s_End	= s_Ini + s_EpochWind;
    
    if s_End > numel(v_Energy)
        s_End   = numel(v_Energy);
    end
    
    v_Perc	= edfcnd(v_Energy(s_Ini:s_End),-inf,[],'method',3);
    v_Val	= v_Perc(:,1);
    v_Perc  = v_Perc(:,2);
    
    s_Index	= find(v_Perc <= s_Percentil,1,'last');
    
    v_Threshold(s_Ini:s_End)= v_Val(s_Index);    
end

v_EnergyThres   = v_Energy >= v_Threshold;


v_WindThres     = [0;v_EnergyThres;0];
v_WindJumps     = diff(v_WindThres);
v_WindJumUp     = find(v_WindJumps==1);
v_WindJumDown   = find(v_WindJumps==-1);
v_WinDist       = v_WindJumDown - v_WindJumUp;

v_WinDistSelect = (v_WinDist > s_MinWind);

v_WindSelect    = find(v_WinDistSelect);

if isempty(v_WindSelect)
    evinfo = [];
    
else
	v_WindIni = v_WindJumUp(v_WindSelect);
	v_WindEnd = v_WindJumDown(v_WindSelect) - 1;
    for k = 1 : length( v_WindIni)
        loc = v_WindIni( k) : v_WindEnd( k);
        evinfo.Location( k, :) = loc( [1, end]);
        [nb_cycle, avg_freq] = EstimateNbcycles( v_FiltSignal( loc), srate);
        evinfo.NoCycles( k, 1) = nb_cycle;
        evinfo.PeakZScore( k, 1) = max( v_Energy( loc));
        evinfo.AvgFreq( k, 1) = avg_freq;
        evinfo.BPFreq( k, :) = cfg.bpfreq;
    end
end



end