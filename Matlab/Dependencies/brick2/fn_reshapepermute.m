function y = fn_reshapepermute(x,varargin)
% function y = fn_reshapepermute(x,dimshift)
% function y = fn_reshapepermute(x,resh,perm[,resh2])
%---
% perform a combination of reshape and permute from a single function
%
% Input:
% - x           input array
% - dimshift    cell array - shows which dimensions will appear together
%               for example, fn_reshapepermute(x,{2 [1 3]}) is equivalent
%               to rehape(permute(x,[2 1 3]),[size(x,2) size(x,1)*size(x,3)])
% - resh, perm, resh2  
%               vectors - reshaping and permuting vectors to be applied one
%               after the other
%
% Output:
% - y           output array

% Thomas Deneux
% Copyright 2010-2012


if nargin==2
    dimshift = varargin{1};
    s = size(x); s(end+1:max([dimshift{:}])) = 1;
    ndimnew = length(dimshift);
    snew = zeros(1,ndimnew);
    for i=1:ndimnew, snew(i) = prod(s(dimshift{i})); end
    y = permute(x,[dimshift{:}]);
    y = reshape(y,snew);
elseif nargin==3 || nargin==4
    [resh perm] = deal(varargin{1:2});
    y = reshape(x,resh);
    y = permute(y,perm);
    if nargin==4
        resh2 = varargin{3};
        y = reshape(y,resh2);
    end
else
    error argument
end