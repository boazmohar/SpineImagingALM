%% fn_readimg

%% Syntax
%  a = fn_readimg(fname[,'permute'])

%% Description
%  read image using imread, and handles additional features:
%  - converts to double
%  - detects if color or gray-scale images (in the last case, use a 2D array per image)
%  - can read a stack of images (returns 3D array)
% 
%  images are read according to x-y convention, use 'permute' flag to use
%  Matlab y-x convention

%% Source
% Thomas Deneux
%
% Copyright 2004-2012
%
