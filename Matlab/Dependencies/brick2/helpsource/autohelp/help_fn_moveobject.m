%% fn_moveobject

%% Syntax
%  dp = fn_moveobject(hobj[,'fast'|('fastcolor',col)][,'latch'][,'point',i])

%% Description
%  moves objects while mouse button is pressed
%  
%  Options
%  - 'fast'      use 'xor' erase mode for faster display updates
%  - 'fastcolor' color to use in 'xor' mode
%  - 'latch'     when button is released, brings objects back to initial
%                position
%  - 'point',i   move only the ith point(s) of line objects
%  - 'twice'     wait for button press+release, or release+pressagain
%  
%  See also fn_buttonmotion

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
