function [filt_vec, amp] = HFOFiltData( vec, srate, param)
%   Band-pass filter and line noise removal using fieldtrip toolbox
% 
%   Input
%       vec,   vector of time series
%       srate, sampling rate in Hz
%       param, data structure with fields:
%               bpfreq,    [lo, hi], band-pass filtering frequency
%               bptype,    'fir' | 'iir', band-pass filtering type
%               notch,     'yes' | 'no', notch filtering
%               notch_freq, 50 | 60 notch frequencies
%               notch_type, 'fir' | 'iir', notch filter type
%               notch_bandwidth, scalar, bandwidth of the notch filter
%
%   Output
%       filt_vec, filtered time series
%            amp, amplitude time series
% 
%   Usage
%       [filt_vec, amp_vec] = HFOFiltData( vec, srate, param);
%
%   Requires fieldtrip toolbox
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
% 

if nargin ~= 3
    error( 'HFOFiltData requires three input arguments');
end

if ~isvector( vec)
    error( 'Vector input required.');
end

if size( vec, 2) == 1
    vec = vec(:)';
end

% Band-pass filtering
if strcmpi( param.bptype, 'fir')
    bptype = 'fir';
else
    bptype = 'but';
end

filt_vec = ft_preproc_bandpassfilter( vec, srate, param.bpfreq, [], bptype);

% Notch filtering
if strcmpi( param.notch, 'yes')  
    if strcmpi( param.notch_type, 'fir')
        notch_type = 'fir';
    else
        notch_type = 'but';
    end

    f = param.notch_freq : param.notch_freq : srate/2 - param.notch_bandwidth; 
    
    for k = 1 : length( f)
        bsfreq = param.notch_freq + 0.5 * [-1, 1] * param.notch_bandwidth;
        filt_vec = ft_preproc_bandstopfilter( filt_vec, srate, bsfreq, [], notch_type);
    end
end

if nargout > 1
    amp = abs( hilbert( filt_vec(:)))';
end


end % function
