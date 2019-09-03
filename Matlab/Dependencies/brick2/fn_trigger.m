function y = fn_trigger(x,indices,varargin)
% function y = fn_trigger(x,indices[,dim][,nframes])
%---
% average according to triggers
%
% Input:
% - x           data to trigger
% - indices     trigger position, in indices coordinates
% - dim         dimension of the data on which to trigger and average [by
%               default, the last dimension is used]
% - nframes     size of the output in this dimension: it should be a
%               2-element vector (how many frames before trigger, and how
%               many frames after trigger), or a scalar (the same number is
%               used both times) [the default is 0, i.e. only averaging at
%               the trigger positions]

% Thomas Deneux
% Copyright 2008-2012

% Input
dim = [];
nframes = 0;
for k=1:length(varargin)
    a = varargin{k};
    if isvector(a) || a>5 || ~isempty(dim)
        nframes = a;
    else
        dim = a;
    end
end
if isempty(dim), dim = ndims(x); end
if isscalar(nframes), nframes = [nframes nframes]; end

% Sizes
s = size(x);
nx = size(x,dim);
ny = sum(nframes)+1;

% Ignore indices too close to edges
bad = (indices<1+nframes(1) | indices>nx-nframes(2));
indices(bad) = [];
nind = length(indices);

% Trigger-averaging
indicesplus = fn_add(indices(:)',(-nframes(1):nframes(2))');
x = reshape(x,[prod(s(1:dim-1)) s(dim:end)]);
y = x(:,indicesplus(:),:);
y = reshape(y,[prod(s(1:dim-1)) ny nind s(dim+1:end)]);
y = mean(y,3);
y = reshape(y,[s(1:dim-1) ny s(dim+1:end)]);

