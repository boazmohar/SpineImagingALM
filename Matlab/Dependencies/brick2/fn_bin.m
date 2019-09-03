function data=fn_bin(data,bins,varargin)
% function data=fn_bin(data,xybin (scalar)[,'same'][,'sum'])
% function data=fn_bin(data,bins (vector)[,'same'][,'sum'])
%---
% bin data according to vector describing which binning to apply for each 
% dimension
% for example, can be used to bin 3D data image (x,y,t) in space and time 
% 'same' flag will result in a binned data of same size as original
% 'sum' flag causes to sum over each mean rather than averaging

% Thomas Deneux
% Copyright 2010-2012

if nargin==0, help fn_bin, return, end

% Input
% (flags)
[dosame dosum] = fn_flags('same','sum',varargin);
% (bin size)
if nargin>2 && ~ischar(varargin{1})
    error('old syntax ''fn_bin(data,xybin,tbin)'' not allowed any more')
end
s = size(data);
nd = ndims(data);
if length(bins)==1 
    if size(data,2)==1
        bins = [bins 1];
    elseif size(data,1)==1
        bins = [1 bins];
    else % special case 'xybin'
        bins=[bins bins];
    end
elseif size(bins,2)==1
    bins = bins';
end
if length(bins)>nd
    nd = length(bins);
    s(end+1:nd) = 1;
elseif length(bins)<nd
    bins(end+1:nd) = 1;
end

% prepare data (resize+reshape)
bins = min(bins,s); % avoid situation where the bin size in some dimension is larger than the data
if dosame
    s2 = ceil(s./bins);
    s1 = s2.*bins;
    if any(s1~=s)
        subs = cell(1,nd);
        for dim=1:nd
            subs{dim} = 1:s(dim);
        end
        data0 = data;
        if numel(data)>1e8, disp('warning: duplicating large array'), end
        data = nan(s1);
        subsasgn(data,substruct('()',subs),data0);
    end
else
    s2 = floor(s./bins);
    s1 = s2.*bins;
    if any(s1~=s)
        subs = cell(1,nd);
        for dim=1:nd
            subs{dim} = 1:s1(dim);
        end
        if numel(data)>1e8, disp('warning: duplicating large array'), end
        data = subsref(data,substruct('()',subs));
    end
end
% (finish if empty data)
if any(s==0), data = zeros(s2); return, end

% bin
stemp = [bins; s2]; stemp = stemp(:)';
data = reshape(data,stemp);
for dim=1:nd
    if bins(dim)>1
        if dosum
            data = sum(data,2*dim-1);
        else
            data = mean(data,2*dim-1);
        end
    end
end

% final reshape+resize
if dosame
    repfact = [bins; ones(1,nd)]; repfact = repfact(:)';
    data = repmat(data,repfact);
    data = reshape(data,s1);
    subs = cell(1,nd);
    for dim=1:nd
        subs{dim} = 1:s(dim);
    end
    data = subsref(data,substruct('()',subs));
else
    data = reshape(data,s2);
end

