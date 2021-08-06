% The automated detectors are modified from RippleLab that
%   Written by:
%   Miguel G. Navarrete Mejia
%   Electrical Engineering MS candidate
%   UNIVERSIDAD DE LOS ANDES
%   Colombia, 2012
%   mnavarretem@gmail.com
% For more details: Navarrete et al., RIPPLELAB: A Comprehensive Application
%       for the Detection, Analysis and Classification of High Frequency Oscillations
%       in Electroencephalographic Signals. Plos One 2016.


function evs = HFOAutoDetect( mat, srate, labels, cfg)
% Automatic HFO events detector
%
% Input
%   mat, channel x sample
%   srate, sampling rate
%   labels, Nx1 cell array of channel labels
%   cfg, configuration data structure
%       ** Required fields:
%       'method', 'ste', 'sll', 'hil' | 'mni', detection methods
%       'bpfreq', [lo, hi] bandpass filter frequency in Hz
%       'bptype', 'fir' | 'iir', filter type
%       'notch', 'yes' | 'no', notch filter
%       'notch_freq', 50 | 60, power line noise frequency in Hz
%       'notch_type', 'fir' | 'iir', notch filter type
%       'notch_bandwidth, scalar, notch filter bandwidth in Hz
%
%       ** Fields dependent on the method
%       If cfg.method == 'hil'
%           'epcoh', scalar, epcoh length in seconds, [] - whole time series
%           'z_thresh', scalar z score threshold
%           'peak_thresh', scalar peak z score threshold, can be the same to cfg.z_thres
%           'nocycle_thresh', set for HFOApp method scalar duration threshold, in seconds
%           'dur_thresh', set for RippleLab method scalar duration threshold, in seconds
%     
%       If cfg.method == 'mni'
%           'epoch', e.g., 10, epoch length in seconds
%           'hf_epcoh,  e.g., 60
%           'hf_percent'  e.g., 95;
%           'min_win', e.g., 0.01;                   % Minimum HFO Time, in seconds   
%           'min_gap', e.g., 0.01;                   % Minimum HFO Gap, in seconds
%           'thresh_percent', 0-100 e.g., 99.9999;         % Threshold precentil
%           'bs_win', e.g., 125;                    % Baseline window , in seconds
%           'bs_shift', e.g. 0.5;                       % Baseline Shift window    
%           'bs_thresh', e.g. 0.67                         % Baseline threshold
%           'bs_min', e.g., 5	                        % Baseline
%     
%       If cfg.method == 'ste'
%           'epoch',  e.g. 180, Cycle Time in seconds
%           'rs_win', e.g.0.003,  scalar, RMS window time in seconds
%           'rs_thresh', e.g.5, Threshold for RMS in standard deviation
%           'rs_min_dur', e.g.0.006,  Min window time for an HFO in seconds
%           'min_dist', e.g.0.010, Min Distance time Betwen two HFO candidates
%           'nocycle_thresh', e.g.6, % Minimum oscillations per interval
%           'peak_thresh', e.g.3, % Threshold for finding peaks 
%     
%       If cfg.method == 'sll'
%           'epoch', e.g. 180, scalar, in seconds
%           'win_len', e.g. 0.005, scalar, length in seconds
%           'percent', e.g. 95.5, 0-100, scalar
%           'min_dur', e.g. 0.012, minimal duration, in seconds
%           'do_equal', e.g. 'yes', 'yes' | 'no' Equalization Filter
% 
%       If cfg.method == other_detector
%           '**'
%
% Output
%   evs, 
%       [], if no HFOs were detected.
%       structure array if HFOs were found
%           evs(x).label = 'A2';
%           evs(x).info.Location, N x 2 location of events
%           evs(x).info.NoCycles, N x 1 vector of number of cycles
%           evs(x).info.PeakZScore, N x 1 vector of peak z score or other measures
%           evs(x).info.AvgFreq, N x 1 average frquency in Hz
%           evs(x).info.BPFreq, N x 2 band-pass frequency [lo, hi] Hz
%       
%     
% Usage
%   evs = HFO_AutoEvents( mat, srate, labels, cfg);
%
% 
%  This file is part of HFOApp.
% 
%       HFOApp is free software: you can redistribute it and/or modify it
%       under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the license, or
%       (at your option) any later version.
% 
%       HFOApp is distributed in the hope that it will be useful, but
%       without any warranty; without even the implied warranty of
%       merchantability or fitness for a particular purpose. See the GNU
%       General Public License for more details.
% 
%       You should have received a copy of the GNU General Public License
%       along with HFOApp. If not, see <http://www.gnu.org/licenses/>.
% 
% G

% validate input


% No events detected
evs = struct([]);

% number of channels and samples
nb_chnls = size( mat, 1);
for chnl_idx = 1 : nb_chnls
    switch cfg.method
        case 'hil'
            % Hilbert detector
            
            % cfg.epcoh = 180;
            % cfg.z_thresh = 5;
            % cfg.peak_thresh = 5;
            % cfg.dur_thresh = 10;
            
            evinfo = HFOAutoDetectorHil( mat( chnl_idx, :), srate, cfg);
            
        case 'mni'
            % MNI detector
            
            % cfg.epoch = 10;
            % cfg.hf_epcoh = 60;
            % cfg.hf_percent  = 95;
            % cfg.min_win = 0.01;
            % cfg.min_gap = 0.01;
            % cfg.thresh_percent = 99.9999;
            % cfg.bs_win = 125;
            % cfg.bs_shift = 0.5;
            % cfg.bs_thresh = 0.67;
            % cfg.bs_min = 5;
            
            evinfo = HFOAutoDetectorMNI( mat( chnl_idx, :), srate, cfg);
            
        case 'ste'
            % Short time energy detector
            
            % cfg.epoch =  180; Cycle Time in seconds
            % cfg.rs_win = 0.003;  scalar, RMS window time in seconds
            % cfg.rs_thresh = 5; Threshold for RMS in standard deviation
            % cfg.rs_min_dur = 0.006;  Min window time for an HFO in seconds
            % cfg.min_dist = 0.010; Min Distance time Betwen two HFO candidates
            % cfg.nocycle_thresh = 6; % Minimum oscillations per interval
            % cfg.peak_thresh = 3; % Threshold for finding peaks
            
            evinfo = HFOAutoDetectorSTE( mat( chnl_idx, :), srate, cfg);
            
        case 'sll'
            % Short line length detector
            
            % cfg.epoch = 180; scalar, in seconds
            % cfg.win_len = 0.005; scalar, length in seconds
            % cfg.percent = 95.5; 0-100, scalar
            % cfg.min_dur = 0.012; minimal duration, in seconds
            % cfg.do_equal = 'yes'; 'yes' | 'no' Equalization Filter
            
            evinfo = HFOAutoDetectorSLL( mat( chnl_idx, :), srate, cfg);
        
        % ****************************************************************
        % case 'your-own-detector'
        % ****************************************************************
        
        
        otherwise
            evinfo = [];
    end
    
    
    if ~isempty( evinfo)
        chnl_ev = [];
        chnl_ev.label = labels{ chnl_idx};
        chnl_ev.Location = evinfo.Location;
        chnl_ev.info = evinfo;
        evs = cat( 1, evs, chnl_ev);
    end
    
end % channel loop

end % main function

