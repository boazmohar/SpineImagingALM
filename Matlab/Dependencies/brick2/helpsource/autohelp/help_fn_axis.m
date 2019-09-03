%% fn_axis

%% Syntax
%  ax = fn_axis([ha,]'tight|image|tightimage',sidefactor,'y0')

%% Description
%  set a nice range to axes
% 
%  Input:
%  - flag            'tight'     stretch the axis as much as possible
%                    'image'     set an equal ratio between x and y
%                    'tightimag' stretch the axis as much as possible, while
%                                maintaining an equal ratio between x and y
%  - sidefactor      scalar or 2-elements vector >1 but close to one: the
%                    difference with 1 indicates the small gap to leave to 
%                    the side
%  - 'y0' flag       set ymin to 0
%  
%  Example:
%    plot(sin(0:.01:100))
%    fn_axis('tight',[1 1.2])

%% Source
% Thomas Deneux
%
% Copyright 2010-2012
%
