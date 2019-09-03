%% fn_fftfrequencies

%% Syntax
%  ff = fn_fftfrequencies(data,fs)
%  [ff1 ff2] = fn_fftfrequencies(data,fs)

%% Description
%  utility to get the vector of frequencies for which fft is computed on
%  data
%  
%  Inputs:
%  - data    vector, or scalar (then it is the number of elements of the
%            data) 
%            data can also be an array, in which case there are 2 outputs
%  - fs      sampling frequency (default=1)

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
