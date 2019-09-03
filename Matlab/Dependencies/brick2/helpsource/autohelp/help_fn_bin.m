%% fn_bin

%% Syntax
%  data=fn_bin(data,xybin (scalar)[,'same'][,'sum'])
%  data=fn_bin(data,bins (vector)[,'same'][,'sum'])

%% Description
%  bin data according to vector describing which binning to apply for each 
%  dimension
%  for example, can be used to bin 3D data image (x,y,t) in space and time 
%  'same' flag will result in a binned data of same size as original
%  'sum' flag causes to sum over each mean rather than averaging

%% Source
% Thomas Deneux
%
% Copyright 2010-2012
%
