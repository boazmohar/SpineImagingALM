%% fn_imageop

%% Syntax
%  a = fn_imageop(a[,par])
%  par = fn_imageop('par')

%% Description
%  Apply a series of transformations to an image
%  
%  Input/Output:
%  - a       2D (or more) array - image (or movie, ...) to be modified
%  - par     structure - parameter of operations to apply; fields are:
%            .bin        binning
%            .xlow       low-pass filter
%            .xhigh      high-pass filter

%% Source
% Thomas Deneux
%
% Copyright 2011-2012
%
