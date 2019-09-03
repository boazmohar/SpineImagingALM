%% fn_savemovie

%% Syntax
%  M = fn_savemovie(a[,fname][,clip][,fps][,zoom][,map][,hf])
%  M = fn_savemovie(a[,'fname',fname][,'clip',clip],...)

%% Description
%  Input:
%  - a       x-y-t array (be aware Matlab convention for images is y-x)
%            or x-y-c-t for true colors (c stands for the 3 color channels,
%            values should be in [0 1] 
%  - fname   file name (movie is saved in file only if specified) [default =
%            30]
%  - clip    a 2-values vector, or 'fit', or '?SD'
%  - fps     frames per second
%  - zoom    zooming value, according to which the movie is either
%            interpolated (zoom>0) or binned (0<zoom<1) or enlarged with
%            "big pixels" (zoom<-1) 
%  - map     nx3 array for the colormap [default = gray(256)]
%  - hf      figure handle (if specified, the movie is played in figure)
% 
%  Note that arguments can be passed in any order, except the first; if
%  there is ambiguity, the function tries to guess which value was entered
%  (for example a scalar value will be assigned to 'fps' if it is >=5, and to
%  'zoom' if it is <5); in order to de-ambiguate, it is possible to preced
%  the value by a flag.
%  e.g. fn_savemovie(rand(3,4,25),'fname','test.avi','zoom',10)
%  
%  Output
%  - M       movie frames (Matlab format)
% 
%  See also fn_readmovie, fn_movie

%% Source
% Thomas Deneux
%
% Copyright 2004-2012
%
