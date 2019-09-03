function b = fn_isemptyc(c)
% function b = fn_isemptyc(c)
%---
% returns an array of logicals of the same size as cell array c indicating
% which elements of c are empty

% Thomas Deneux
% Copyright 2011-2012

b = false(size(c));
for k=1:numel(c), b(k) = isempty(c{k}); end
