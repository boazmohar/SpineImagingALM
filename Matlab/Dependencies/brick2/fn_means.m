function y = fn_means(varargin)
% function y = fn_means(x1,x2,x3...)
%---
% average all arguments (they need to be same size)

% Thomas Deneux
% Copyright 2004-2012

if nargin<1, help fn_means, return, end

doconv = ~fn_ismemberstr(class(varargin{1}),{'double' 'single'});

n = 0;
y = [];
for i=1:nargin
    if ~isempty(varargin{i})
        n = n+1;
        xi = varargin{i};
        if doconv, xi = single(xi); end
        if isempty(y)
            y = xi;
        else
            y = y+xi;
        end
    end
end
y = y/n;