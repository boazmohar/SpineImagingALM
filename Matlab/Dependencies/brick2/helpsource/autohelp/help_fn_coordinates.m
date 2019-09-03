%% fn_coordinates

%% Syntax
%  [orig scale] = fn_coordinates([handle,]refchangetag)
%  x2 = fn_coordinates([handle,]refchangetag,x1,'position')
%  v2 = fn_coordinates([handle,]refchangetag,v1,'vector')

%% Description
%  Change coordinates between screen (pixels), figure, axes
%  The "referential change tag" will be 's2f', 'f2a', ...
%  Possible referentials are:
%  - s       screen, bottom-left origin, pixel units
%  - f       figure, bottom-left origin, normalized units
%  - g       figure, bottom-left origin, pixel units
%  - a       axes, axes units
%  - b       axes, bottom-left origin, pixel units
%  - c       axes, bottom-left origin, normalized units
%  
%  Ex: the axes coordinates of the point at position 10 pixels to the right 
%  and the top of the axes origin in figure 1 are
%  xy = fn_coordinates(1,'b2a',[10 10],'position')

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
