%% fn_framedisplay

%% Syntax
% [y =]  fn_framedisplay(x[,clip][,'singleaxes|multaxes'][,hf|ha][,'display'][,ncol])

%% Description
%  Input (except for x, order of inputs can be changed):
%  - x       3D or 4D data; NOT using Matlab image convention (i.e. first
%            dim is x and second dim is y) 
%  - clip    clipping range
%  - axflag  'singleaxes' will create a large image by putting each frame
%            next to each other and display this image
%            'multaxes' will display each frame in a separate axes
%  - hf|ha   handle of figure (implicit 'multaxes' option) or axes (implicit
%            'singleaxes' option)
%  - 'display'   by default, nothing is displayed if an output is requested,
%            this option is used to display even if there is an output
%  - ncol    specify the number of columns to use
% 
%  Output:
%  - y       a large image obtained by concatenating all images, note that
%            compared to the large image created with the 'singleaxes'
%            option, this image is slighlty larger since there is a 1-pixel
%            wide separation between all frames

%% Source
% Thomas Deneux
%
% Copyright 2011-2012
%
