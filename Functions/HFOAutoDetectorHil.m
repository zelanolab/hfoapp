function evinfo = HFOAutoDetectorHil( v, srate, cfg)
% Hilbert detector of high frequency oscillation events
% 
% Input
%   v, row vector
%   srate, sampling rate in Hz
%   cfg, data structure
%       'epcoh', scalar, epcoh length in seconds, [] - whole time series
%       'z_thresh', scalar z score onset threshold
%       'peak_thresh', scalar peak z score threshold, can be the same to cfg.z_thresh
%       'dur_thresh', scalar duration threshold, in seconds
%       'nocycle_thresh', [] (default), scalar, set this to filter events
%       'avg_method', 'mov' | 'epoch'
% 
% Output
%   evinfo, [] if no events else
%       data structure with the following fields
%       
%       'Location', N x 2 matrix, event location
%       'NoCycles', N x 1 vector, number of cycles
%       'PeakZScore', N x 1 vector, peak z score
%       'AvgFreq', N x 1 vector, average frequency in Hz
%       'BPFreq', N x 2 matrix, [lo, hi] frequency boundary in Hz
% 
% G


if nargin < 3
    error( 'Not engouh input arguments.');
end


nb_pnts = length( v);

if ~isempty( cfg.epoch)
    epoch_len = round( cfg.epoch * srate);
    epoch_idx = buffer( 1 : nb_pnts, epoch_len, 0, 'nodelay');
    nb_epochs = size( epoch_idx, 2);
    
else
    nb_epochs = 1;
    epoch_idx = (1 : nb_pnts)';
    epoch_len = nb_pnts;
end

% Compute envelope of band-pass filtered time series
[filt_v, amp] = HFOFiltData( v, srate, cfg);
if strcmpi( cfg.avg_method, 'mov')
    zz = (amp - movmean( amp, epoch_len)) ./ movstd( amp, epoch_len);
end

evinfo = [];
evinfo.Location = []; % N x 2, event location
evinfo.NoCycles = []; % N x 1, number of cycles
evinfo.PeakZScore = []; % N x 1, peak z score
evinfo.AvgFreq = []; % N x 1, average frequency
evinfo.BPFreq = []; % N x 2, [lo, hi] frequency boundary

% Loop through epochs
cnt = 0;
for x = 1 : nb_epochs
    % z-score normalization
    loc = epoch_idx( epoch_idx( :, x) > 0, x);
    
    if strcmpi( cfg.avg_method, 'epoch')
        z = zscore( amp( loc));
    else
        z = zz( loc);
    end
    
    % threshold to find events
    bin_z = z >= cfg.z_thresh;
    if ~any( bin_z)
        continue;
    end
    
    ind = find( bin_z);
    ind = ind(:);    
    b = find( diff( ind) > 1);
    b = b(:);
    on_loc = ind( [1; b+1]);
    off_loc = ind( [b; length( ind)]);
    
    % maximal z value within each event
    c = cellfun( @(s, e) {s : e}, num2cell( on_loc), num2cell( off_loc));
    peakz = cellfun( @(y) max( z(y)), c, 'uniformoutput', false);
    peakz = cell2mat( peakz);
    
    % duration
    dur = (off_loc - on_loc) / srate;
    
    if isfield( cfg, 'dur_thresh') && ~isempty( cfg.dur_thresh)
        ev_idx = find( dur >= cfg.dur_thresh & peakz >= cfg.peak_thresh);
    else
        ev_idx = find( peakz >= cfg.peak_thresh);
    end
    
    for k = 1 : length( ev_idx)
        loc = [on_loc( ev_idx( k)), off_loc( ev_idx( k))] + epoch_idx( 1, x) - 1;
        
        % Estimate number of cycles and average frequency
        [nb_cycles, avg_freq] = EstimateNbcycles( filt_v(loc(1) : loc(2)), srate);
        
        if isfield( cfg, 'nocycle_thresh') && ~isempty( cfg.nocycle_thresh)
            if nb_cycles < cfg.nocycle_thresh
                continue;
            end
        end
            
        cnt = cnt + 1;        
        evinfo.Location = cat( 1, evinfo.Location, loc);
        evinfo.NoCycles = cat( 1, evinfo.NoCycles, nb_cycles);
        evinfo.PeakZScore = cat( 1, evinfo.PeakZScore, peakz( ev_idx( k)));
        evinfo.AvgFreq = cat( 1, evinfo.AvgFreq, avg_freq);
        evinfo.BPFreq = cat( 1, evinfo.BPFreq, cfg.bpfreq);
    end    
end % epoch loop

if cnt < 1
    evinfo = [];
end

end % hilbert detector
