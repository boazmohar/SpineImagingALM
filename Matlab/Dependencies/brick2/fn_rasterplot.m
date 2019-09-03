function fn_rasterplot(varargin)
% function fn_rasterplot(t[,i,...])
% function fn_rasterplot(c,...)
%---
% Input:
% - t       time points
% - i       index of raster line they belong to
% or:
% - c       cell array with time points for each separate line (cell array
%           elements must be column vectors)

% Thomas Deneux
% Copyright 2008-2012

if nargin==0, help tps_rasterplot, return, end
if iscell(varargin{1})
    t = varargin{1}; varargin(1)=[];
    i = t; for k=1:numel(i), i{k}(:)=k; end
    try
        t = cat(1,t{:});
        i = cat(1,i{:});
    catch %#ok<CTCH>
        t = cat(2,t{:})';
        i = cat(2,i{:})';
    end
else
    t = varargin{1}; varargin(1)=[];
    if nargin<2
        i = ones(1,length(t));
    else
        i = varargin{1}; varargin(1)=[];
        if isscalar(i), i = ones(1,length(t))*i; end
    end
end

t = repmat(t(:)',3,1);
i = fn_add([0 .5 NaN]',i(:)');

plot(t(:),i(:),varargin{:})