function [nb_cycle, avg_freq] = EstimateNbcycles( v, srate)
% Estimate average frequency and srate based on peaks.
% 
% Input
%   v, vector
%   srate, sampling rate in Hz
% 
%   minimal peak distance is set to 2/srate
% 
% Output
%   nb_cycle, number of cycles
%   avg_freq, average frequency: length / avgerage_cycle
% 
% G

warning('off', 'signal:findpeaks:largeMinPeakHeight');

nb_cycle = 0;
avg_freq = 0;

thresh = 0;
if ~any( v > thresh) || length( v) < 3
    return;
end

try
    [~, pklocs] = findpeaks( v, 'MinPeakHeight', thresh, 'MinPeakDistance', 1/srate);
    
    if length( pklocs) > 1
        cycle_dur = mean( diff( pklocs));
        nb_cycle = length( v) / cycle_dur;
        avg_freq = srate / cycle_dur;
    end
end

end % function