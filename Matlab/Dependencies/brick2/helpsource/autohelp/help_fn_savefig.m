%% fn_savefig

%% Syntax
%  fn_savefig([hf][,fnames,...[,methods...]][,scaling])

%% Description
% 
%  Input:
%  - hf          figure handles
%  - fnames      file name where to save the figure
%  - methods     for example: 'jpg', 'eps', 'fig', 'png'
%  - scaling     scalar: how much to increase the size of the figure (use
%                small values to increase font size)
% 
%  For example:
%  * fn_savefig(1,'myfigure','png','eps') saves figure 1 in files
%  'myfigure.png' and 'myfigure.eps' 
%  * fn_savefig(1:3,'a.png','b.jpg','c.fig',2) save figure 1 in 'a.png',
%  figure 2 in 'b.jpg' and figure 2 in 'c.fig'; image sizes are double as
%  default
%  * fn_savefig(1) asks you where to save figure 1

%% Source
% Thomas Deneux
%
% Copyright 2003-2012
%
