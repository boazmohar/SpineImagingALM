%% fn_reshapepermute

%% Syntax
%  y = fn_reshapepermute(x,dimshift)
%  y = fn_reshapepermute(x,resh,perm[,resh2])
% perform a combination of reshape and permute from a single 

%% Description
% 
%  Input:
%  - x           input array
%  - dimshift    cell array - shows which dimensions will appear together
%                for example, fn_reshapepermute(x,{2 [1 3]}) is equivalent
%                to rehape(permute(x,[2 1 3]),[size(x,2) size(x,1)*size(x,3)])
%  - resh, perm, resh2  
%                vectors - reshaping and permuting vectors to be applied one
%                after the other
% 
%  Output:
%  - y           output array

%% Source
% Thomas Deneux
%
% Copyright 2010-2012
%
