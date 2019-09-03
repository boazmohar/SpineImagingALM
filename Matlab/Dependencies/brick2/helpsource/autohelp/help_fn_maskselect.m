%% fn_maskselect

%% Syntax
%  mask = fn_maskselect(image[,mouseflag][,dorepeat])

%% Description
%  
%  Input:
%  - image       2D array
%  - mouseflag   'rect', 'poly' [default], 'free', 'ellipse'
%  - dorepeat    select multiple regions? [default = false]
%  
%  Output:
%  - mask        logical array the same size of image indicating interior of
%                the mask

%% Source
% Thomas Deneux
%
% Copyright 2011-2012
%
