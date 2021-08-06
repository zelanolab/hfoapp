%***************************************************************************
% Modified from f_GaborTransformWait.m in RippleLab
%
% f_GaborAWTransformMatlab Description:
% This function calculates the Wavelet Transform using a Gaussian modulated
% window (Gabor Wavelet).
% The Sample Rate of the input signal is considered in order to compute
% the transform.
% Author: Mario Valderrama
%
% ps_StDevCycles: In the wavelet transform, the scale corresponding to each
% frequency (in v_FreqAxis) defines the value (in seconds) of the
% gaussian's standard deviation used for the calculation of the transform
% at every one of these frequencies. This standard deviation value must be
% big enough to cover at least one or more complete cycles (or periods) of
% the oscillation at each considered frequency. Thus, this parameter
% defines the number of cycles you want to include for the transform at
% every frequency.
%
% ps_Magnitudes: Set to 1 (default) if the magnitudes of the coefficients
% must be returned; 0 for analytic values (complex values).
%
% ps_SquaredMag: Set to 1 if the magnitudes of the coefficients divided by
% the squared of the corresponding scale must by power to 2
%
% ps_MakeBandAve: Set to 1 if instead of returning a matrix with all values
% in the time-frequency map, the function returns just a vector with the
% average along all the frequency scales for each time moment.
%
% ps_TimeStep: Time step between values that are going to be kept in the
% output matrix. Each time moment is the average of the previous values
% according to the size of the window defined by this parameter.
%
% Outputs:
% m_GaborWT: Matrix containing the scalogram. Time in rows, frequency in
% colums. Frequencies in descending order
% v_TimeAxis: Array containing the time axis values (second units)
% v_FreqAxis: Array containing the frequency axis values in descending
% order (Hz units)
%
% Miguel Navarrete, Catalina Alvarado-Rojas,  Michel Le Van Quyen, Mario Valderrama.
%  RIPPLELAB: A Comprehensive Application for the Detection,
%  Analysis and Classification of High Frequency Oscillations in Electroencephalographic Signals
%  Plos One 2016,  https://doi.org/10.1371/journal.pone.0158276
%***************************************************************************




function [m_GaborWT, t, f] = HFOSpectrogram( pv_Signal, srate, varargin)
% Time-frequency analysis
%
% Input
%   pv_Signal, vector of time series
%   srate, sampling rate in Hz
%   Key-value pairs
%       'foi', vector of frequency of interests
%       'toi', vector of time of interest
%       'ps_StDevCycles', see f_GaborTransformWait.m
%       'ps_Magnitudes', see f_GaborTransformWait.m
%       'ps_SquaredMag', see f_GaborTransformWait.m
%       'ps_TimeStep', see f_GaborAWTransformMatlab.m
%
% Output
%   See f_GaborAWTransformMatlab.m
% 
%  This file is part of HFOApp.
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

options = struct( 'foi', [], ... % foi = linspace( 30, 250, 100);
    'toi', [],... % vector of time of interest
    'ps_StDevCycles', 3,...
    'ps_Magnitudes', 1,...
    'ps_SquaredMag', 0,...
    'ps_TimeStep', []);

if nargin < 2 || ~isvector( pv_Signal)
    error( 'Not enough input arguments');
end

option_names = fieldnames( options);
nbargus = length( varargin);
if mod( nbargus, 2) ~= 0
    error('Key-Value pairs needed.')
end

for pair = reshape( varargin, 2, [])
    in_name = lower( pair{1});
    if any( strcmpi( in_name, option_names))
        options.( in_name) = pair{2};
    else
        error('%s is not a recognized parameter name.', in_name)
    end
end

ps_Magnitudes = options.ps_Magnitudes;
ps_SquaredMag = options.ps_SquaredMag;
ps_StDevCycles = options.ps_StDevCycles;
ps_TimeStep = options.ps_TimeStep;

pv_Signal = pv_Signal(:);
options.foi = transpose( options.foi(:));
f = fliplr( options.foi);

    %if mod(numel(pv_Signal), 2) == 0
    %    pv_Signal = pv_Signal(1:end - 1);
    %end

nbsamples = length( pv_Signal);
v_WAxis = linspace( 0, 2*pi, nbsamples + 1);
v_WAxis = v_WAxis.* srate;
s_HalfLen = floor(nbsamples / 2) + 1;
v_WAxisHalf = v_WAxis( 1 : s_HalfLen);

t = linspace( 0, (nbsamples - 1)/srate, nbsamples);
s_SampAve = max( [1, round(ps_TimeStep * srate)]);
v_SampAveFilt = [];
if s_SampAve > 1
    v_IndSamp = 1 : s_SampAve : nbsamples;
    t = t( v_IndSamp);
    v_SampAveFilt = ones( s_SampAve, 1);
end

v_fft = fft( pv_Signal, nbsamples);
nbfreqs = length( options.foi);
m_GaborWT = zeros( nbfreqs, nbsamples);
s_FreqInd = 0;
for s_FreqCounter = f
    
    v_WinFFT = zeros( nbsamples, 1);
    s_StDevSec = (1 / s_FreqCounter) * ps_StDevCycles;    
    v_WinFFT( 1 : s_HalfLen) = exp(-0.5.* ...
        realpow( v_WAxisHalf - (2.* pi.* s_FreqCounter), 2) .* (s_StDevSec.^ 2));
    v_WinFFT = v_WinFFT .* sqrt( nbsamples) ./ norm(v_WinFFT, 2);
    
    s_FreqInd = s_FreqInd + 1;
    if s_SampAve > 1
        v_GaborTemp = zeros( nbsamples + (s_SampAve - 1), 1);
        v_GaborTemp( s_SampAve : end) =  ifft( v_fft.* v_WinFFT) ./ sqrt(s_StDevSec);
        
        if ps_Magnitudes
            v_GaborTemp = abs( v_GaborTemp);
        end
        
        if ps_SquaredMag
            v_GaborTemp = v_GaborTemp.^2;
        end
        
        v_GaborTemp( 1 : (s_SampAve - 1)) = ...
            flipud(v_GaborTemp(s_SampAve + 1:2 * s_SampAve - 1));
        
        v_GaborTemp = filter( v_SampAveFilt, 1, v_GaborTemp)./ s_SampAve;
        v_GaborTemp = v_GaborTemp( s_SampAve:end);
        m_GaborWT(s_FreqInd, :) = v_GaborTemp( v_IndSamp);
        
    else
        m_GaborWT( s_FreqInd, :) = ifft( v_fft .*  v_WinFFT) ./ sqrt(s_StDevSec);
    end
end

m_GaborWT = flipud( m_GaborWT);
f = fliplr( f);

end % function