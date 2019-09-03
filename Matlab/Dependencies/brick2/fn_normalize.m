function x = fn_normalize(x,dim,flag)
% function x = fn_normalize(x,dim,flag)
%---
% Input:
% - x       array
% - dim     dimensions on which to operate (can be a cell array for
%           multiple actions)
% - action  'div'[default], 'sub', 'std' (can be a cell array for
%           multiple actions), 'detrend'

% Thomas Deneux
% Copyright 2005-2012

if nargin<2, dim=1; end
if nargin<3, flag='div'; end

if isempty(x), return, end

if iscell(dim)
    nrep = length(dim);
    if ~iscell(flag), [tmp{1:nrep}] = deal(flag); flag=tmp; end
    for i=1:nrep
        x = fn_normalize(x,dim{i},flag{i});
    end
    return
end
    
nd = length(dim);
if ~strcmp(class(x),'single'), x = double(x); end
m = x;
for k=1:nd
    m = mean(m,dim(k));
end
switch flag
    case {'div','/'}
        x = fn_mult(x,1./m);
    case {'sub','-'}
        x = fn_add(x,-m);
    case 'std'
        x = fn_add(x,-m);
        m2 = x.^2;
        for k=1:nd
            m2 = mean(m2,dim(k));
        end
        x = fn_mult(x,1./sqrt(m2));
    case 'detrend'
        if nd~=1, error('detrend should be along only one dimension'), end
        nt = size(x,dim);
        tt = shiftdim((1:nt)',1-dim); tt = tt-mean(tt); tt = tt/sqrt(sum(tt.^2));
        a = sum(fn_mult(x,tt),dim);
        x = x - fn_add(fn_mult(a,tt),m);
    otherwise
        error(['unknown flag ''' flag ''''])
end

