%% fn_errorbar

%% Syntax
%  hl = fn_errorbar([x,]y,e,line options[,flag])
%  hl = fn_errorbar([x,]y,line options[,flag])

%% Description
%  ym and ystd can be vectors or arrays
%  if e is not supplied, y is redefined and e is defined as the mean and
%  std/sqrt(n) of y along the 2d dimension (if y is a 2D array) or the 3d
%  dimension (if y is a 3D array)
% 
%  flag can be 'bar' or 'patch'
% 
%  uses bar display instead of plot if 'bar' flag is specified

%% Source
% Thomas Deneux
%
% Copyright 2006-2012
%
