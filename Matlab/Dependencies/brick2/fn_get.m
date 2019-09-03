function x = fn_get(h,f,flag)
% function x = fn_get(h,f[,'cell|struct'])
%--- 
% returns values for properties f of objects with handle h
% output is a cell array or a structure according to flag
% if 'flag' is not specified, output is a structure if f is a cell array of
% strings
% 
% See also fn_set

% Thomas Deneux
% Copyright 2007-2012

if nargin==0, help fn_get, return, end

% Input
if nargin<3
    if iscell(f), flag='struct'; else flag='cell'; end
end
if ~iscell(f)
    f = cellstr(f);
end
nobj   = length(h);
nfield = length(f);

% Output
x = cell(nobj,nfield);
for i=1:nobj
    for j=1:nfield
        x{i,j} = get(h(i),f{j});
    end
end
if strcmp(flag,'struct')
    x = cell2struct(x,f,2);
end
