function [data, status] = HFOLoadData( )
% Load .trc or .mat when the Open Menu is pressed in HFHFOApp
% 
% Input
%   None.
%   
% Supported format
%   Matlab file '*.mat'
%   European Data Format (*.edf)
%   micromed format (*.trc)
%
%
% Output
%   data, a structure with following fields
%               mat,  channel x samples matrix
%             srate,  scalar sampling rate in Hz
%            labels,  Nx1 cell array of channel labels
%       start_clock,  clock time of start, Matlab datetime, 'yyyy-MM-dd HH:mm:ss'
%              file,  full-path-to-file
%
%   status, status of file loading
%       1, succeed
%       0, unsucceeded
%
%   [data, status] = HFOLoadData( );
% 
% 
%  This file is part of HFOApp.
%       HFHFOApp is free software: you can redistribute it and/or modify it
%       under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the license, or
%       (at your option) any later version.
% 
%       HFHFOApp is distributed in the hope that it will be useful, but
%       without any warranty; without even the implied warranty of
%       merchantability or fitness for a particular purpose. See the GNU
%       General Public License for more details.
% 
%       You should have received a copy of the GNU General Public License
%       along with HFHFOApp. If not, see <http://www.gnu.org/licenses/>.
% 
% 

data = [];

[fname, p] = uigetfile( {'*.*'}, 'Please select file to load.');
fname_in_full = fullfile( p, fname);

try
    % *trc files
    [~, ~, ext] = fileparts( fname_in_full);
    
    switch lower( ext)
        case '.trc' %% TRC file
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Copied from fieldtrip toolbox
            % mat = ft_read_data( fname_in_full);
            % hdr = ft_read_header( fname_in_full);            
            % fopen_or_error.m read_micromed_trc.m
            orig = read_micromed_trc( fname_in_full);
            
            hdr             = [];
            hdr.Fs          = orig.Rate_Min; % FIXME is this correct?
            hdr.nChans      = orig.Num_Chan;
            hdr.nSamples    = orig.Num_Samples;
            hdr.nSamplesPre = 0; % continuous
            hdr.nTrials     = 1; % continuous
            hdr.label       = cell(1,hdr.nChans);
            % give this warning only once
            hdr.label  = {orig.elec.Name};
            hdr.chanunit = {orig.elec.Unit};
            hdr.subjectname = orig.name;
            %warning('using a modified read_micromed_trc() function');
            
            % this should be a column vector
            hdr.label = hdr.label(:);
            % remember the original header details
            hdr.orig = orig;
            
            hdr.Fs          = double(hdr.Fs);
            hdr.nSamples    = double(hdr.nSamples);
            hdr.nSamplesPre = double(hdr.nSamplesPre);
            hdr.nTrials     = double(hdr.nTrials);
            hdr.nChans      = double(hdr.nChans);
            
            data.mat = read_micromed_trc( fname_in_full, 1, hdr.nSamples);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % To get starting clock time
            sFile = in_fopen_micromed( fname_in_full);
            t = sFile.header.acquisition;
            start_time = [num2str( t.year), '-', num2str( t.month), '-', num2str( t.day),...
                ' ', num2str( t.hour), ':', num2str( t.min), ':', num2str( t.sec)];
            
            % start_time = sprintf( '%d-%d-%d %d:%d:%d', t.year, t.month, t.day, t.hour, t.min, t.sec);
            start_clock = datetime( start_time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss' );
            start_clock.Format = 'yyyy-MM-dd HH:mm:ss';
            
            data.labels = cellfun( @(x) deblank(x), hdr.label, 'uniformoutput', false);
            data.srate = hdr.Fs;
            data.start_clock = start_clock;
            data.file = fname_in_full;
            status = 1;            
            
        case '.mat' %% Matlab file            
            data = load( fname_in_full);
            data.file = fname_in_full;
            status = 1;            
            
        case '.edf'   
            hdr = read_edf( fname_in_full);
            data.mat = read_edf( fname_in_full, hdr, 1, hdr.nSamples, 1:hdr.nChans);
  
            data.labels = cellfun( @(x) deblank(x), hdr.label, 'uniformoutput', false);
            data.srate = hdr.Fs;
            
            v = hdr.orig.T0;
            start_time = sprintf( '%d-%d-%d %d:%d:%d', v(1), v(2), v(3), v(4), v(5), v(6));            
            start_clock = datetime( start_time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss' );
            data.start_clock = start_clock;
            data.file = fname_in_full;
            
            status = 1;            

        otherwise %% Add more supported file type here
            status = 0;
    end
    
    
catch
    data = [];
    status = 0;
end


end % function