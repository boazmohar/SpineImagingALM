%% fn_filter

%% Syntax
%  y=fn_filter(X0,y[,highpassflag])

%% Description
%  remove confound subspace defined by the matrix X0 to y :
%  y = y - X0*X0'*y
%  if y is a vector, returns column or row vector according to input
%  if y is a matrix, filter is applied to columns
%  if highpassflag is set to false [default=true], only the confounds are kept

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
