function [colors ncol] = fn_colorset(k)
% function [colors ncol] = fn_colorset
% function color = fn_colorset(k)
%---
% set of colors that is larger than the Matlab default 'ColorOrder'
% property of axes
%
% if k is specified, returns the kth color, otherwise returns the set of
% colors and the length of this set

% Thomas Deneux
% Copyright 2008-2012

colors = [0 0 1 ; 0 .75 0 ; 1 0 0 ; ...
    0 .75 .75 ; .75 0 .75 ; .75 .75 0 ; 0 0 0 ; ...
    .75 .35 0 ; 0 1 0 ; 0 .3 0 ; .3 0 0 ; .3 0 .5];
ncol = size(colors,1);
if nargin>=1, colors = colors(1+mod(k-1,ncol),:); end
