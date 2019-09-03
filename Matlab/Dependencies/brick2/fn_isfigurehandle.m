function [b c] = fn_isfigurehandle(x)
% function [isfig isnewfig] = fn_isfigurehandle(x)

% Thomas Deneux
% Copyright 2008-2012

if ~isscalar(x)
    b = false(size(x));
    for i=1:numel(x), b(i) = fn_isfigurehandle(x(i)); end
    return
end

b = (ishandle(x) && strcmp(get(x,'type'),'figure'));
c = false;
if ~b && isnumeric(x) && mod(x,1)==0 && x>0
    b = true;
    c = true;
    figure(x)
end