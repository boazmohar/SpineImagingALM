function c = fn_ismemberstr(a,b)
% function c = fn_ismemberstr(a,b)
%---
% same as ismember(a,b), but much faster!!, for cell arrays of strings

% Thomas Deneux
% Copyright 2007-2012

if ~iscell(a)
    a = {a}; 
    if ~iscell(b)
        error('one of the two argument must be a cell array of strings, use ismember instead')
    end
elseif ~iscell(b)
    b = {b}; 
end

c = false(size(a));
for i=1:numel(a)
    ai = a(i);
    for j=1:numel(b)
        if strcmp(ai,b{j}), c(i)=true; break, end
    end
end

