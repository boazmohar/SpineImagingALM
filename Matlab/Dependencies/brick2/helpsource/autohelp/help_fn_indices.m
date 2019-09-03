%% fn_indices

%% Syntax
%  globi = fn_indices(sizes|array,i,j,k,...)
%  globi = fn_indices(sizes|array,ijk)
%  [i j k...] = fn_indices(sizes|array,globi)
%  ijk = fn_indices(sizes|array,globi)

%% Description
%  converts between global and per-coordinate indices
%  
%  Input:
%  - array   array - use conversion for array of this size (we note nd its
%            number of dimensions)
%  - size    vector - sizes of the array
% 
%  Input/Output:
%  - i,j,k   scalar or vectors of the same length N - per-coordinates indices
%  - ijk     nd vector or nd x N array - per-coordinates indices
%  - globi   scalar or vector of length N - global indices
% 
%  See also fn_imvect

%% Source
% Thomas Deneux
%
% Copyright 2004-2012
%
