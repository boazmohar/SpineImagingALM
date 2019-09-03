%% fn_sensor

%% Syntax
%  S = fn_sensor(propname1,propvalue1,...)

%% Description
%  creates a 'sensor' control, object optimally designed for controlling
%  the clipping range of an image
%  to manually change the clipping value, click on the control and move
%  the pointer while keeping the button pressed; there are 3 different
%  ways of controling the change depending on which button was pressed:
%  - left button     move horizontally to change contrast, and
%                    vertically to change luminosity
%  - middle button   move horizontally to change the minimum of the
%                    clipping range, and vertically to change the
%                    maximum
%  - right button    changes the contrast, but not the mean value of the
%                    clipping range

%% Source
% Thomas Deneux
%
% Copyright 2008-2012
%
