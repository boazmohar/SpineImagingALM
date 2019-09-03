function varargout = fn_fftfrequencies(data,fs)
% function ff = fn_fftfrequencies(data,fs)
% function [ff1 ff2] = fn_fftfrequencies(data,fs)
%---
% utility to get the vector of frequencies for which fft is computed on
% data
% 
% Inputs:
% - data    vector, or scalar (then it is the number of elements of the
%           data) 
%           data can also be an array, in which case there are 2 outputs
% - fs      sampling frequency (default=1)

% Thomas Deneux
% Copyright 2007-2012

if nargin<1, help fn_fftfrequencies, return, end
if nargin<2, fs=1; end

if isscalar(data)
    n = data;
elseif isvector(data)
    n = length(data);
else
    n = size(data);
    if isscalar(fs), fs=[fs fs]; end
end

if length(n)==1
    ff = fs * (0:(n-1))/n;
    varargout = {ff};
else
    ff1 = fs(1) * (0:(n(1)-1))/n(1);
    ff2 = fs(2) * (0:(n(2)-1))/n(2);
    if nargout==0
        varargout = {{ff1 ff2}};
    else
        varargout = {ff1 ff2};
    end
end