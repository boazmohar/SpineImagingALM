%% fn_mouse

%% Syntax
%  poly = fn_mouse([axes handle],'point|cross|poly|free|ellipse'[,msg])
%  [x y] = fn_mouse([axes handle],'point|cross'[,msg])
%  rect = fn_mouse([axes handle],'rect'[,msg])
%  [center axis e] = fn_mouse([axes handle],'ellipse'[,msg])
% multi-s  using mouse events

%% Description
%  mode defines action:
%  'point'       [default] get coordinates on mouse click
%  'cross'       get coordinates on mouse click - use cross pointer
%  'rect'        get a rectangle selection (format [xstart ystart xsize ysize])
%  'rectangle'   get a rectangle selection (format [x1 x2 x3 x4; y1 y2 y3 y4])
%  'poly'        polygone selection
%  'free'
%  'ellipse'
%  options: (ex: 'rect+', 'poly-@:.25:')
%  +     selection is drawn (all modes)
%  -     use as first point the current point in axes (rect, poly, free, ellipse)
%  @     plots open line instead of closed polygon (poly, free, ellipse)
%  :num: interpolates line with one point every 'num' pixel (poly, free, ellipse)
%        for 'ellipse' mode, if :num; is not specified, output is a cell
%        array {center axis e} event in the case of only one outpout argument
%  
%  See also fn_maskselect

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
