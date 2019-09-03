function y = fn_interleave(dim,varargin)
% function y = fn_interleave(dim,x1,x2,...)
%---
% Similar syntax to Matlab function CAT, but interleaves the data

% Thomas Deneux
% Copyright 2012-2012

% Input
x  = varargin;
nx = length(x);
siz  = size(x{1});

% Reshape 
siz2 = [siz(1:dim-1) 1 siz(dim:end)];
for i=1:nx
    xi = x{i};
    if any(size(xi)~=siz), error('size mismatch'), end
    x{i} = reshape(xi,siz2);
end

% Concatenate
y = cat(dim,x{:});

% Reshape
siz3 = [siz(1:dim-1) siz(dim)*nx siz(dim+1:end)];
y = reshape(y,siz3);

