function siz = fn_pixelsize(hobj)
% function siz = fn_pixelsize(hobj)
%---
% returns the width and height in pixels of any object without needing to
% change any units values

% Thomas Deneux
% Copyright 2011-2012

pos = fn_pixelpos(hobj);
siz = pos(3:4);