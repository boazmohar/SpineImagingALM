function fn_lines(a,b,varargin)
% function fn_lines('x|y',coordinates[,ha][,line options])
% function fn_lines(xcoordinates,ycoordinates[,ha][,line options])
%---
% draw a series of vertical and/or horizontal lines

% Thomas Deneux
% Copyright 2004-2012

% Input
dox = true; doy = true;
if ischar(a)
    switch a
        case 'x'
            x = b;
            doy = false;
        case 'y'
            dox = false;
            y = b;
    otherwise
        error 'type must be ''x'' or ''y'''
    end
else
    x = a;
    y = b;
end
if ishandle(varargin{1})
    ha = varargin{1};
    varargin(1) = [];
else
    ha = gca;
end

% Go
ax = axis(ha);
if dox
    for k=1:length(x)
        line([1 1]*x(k),ax([3 4]),varargin{:})
    end
end
if doy
    for k=1:length(y)
        line(ax([1 2]),[1 1]*y(k),varargin{:})
    end
end


