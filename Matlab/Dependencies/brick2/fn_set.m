function fn_set(h,f,x)
% function fn_set(h,f,x)
% function fn_set(h,x)
%--- 
% performs set(h(i),f{j},x{i,j}) for every possible i and j
% if x is a structure, performs set(h(i),f{j},x(i).(f{j})), and f input is
% facultative
% 
% See also fn_get

% Thomas Deneux
% Copyright 2007-2012

if nargin==0, help fn_set, return, end

% Input
if nargin==2
    x = f;
    f = fieldnames(x);
elseif ~iscell(f)
    f = cellstr(f);
end
nobj   = length(h);
nfield = length(f);
if ~iscell(x) && ~isstruct(x)
    x = {x};
end
structflag = isstruct(x);
if structflag
    if isscalar(x)
        x = repmat(x,nobj,1);
    end
elseif isvector(x)
    if isscalar(x)
        x = repmat(x,nobj,nfield);
    elseif length(x)==nobj
        x = repmat(x(:),1,nfield);
    end
end
if fn_switch(structflag,length(x)~=nobj,size(x)~=[nobj nfield])
    error('dimension mismatch')
end

% Output
for i=1:nobj
    for j=1:nfield
        fj = f{j};
        if structflag, xij = x(i).(fj); else xij = x{i,j}; end  
        if isprop(h(i),fj), set(h(i),fj,xij); end
    end
end
