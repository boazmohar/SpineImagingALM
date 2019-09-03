%% fn_sliderenhance

%% Syntax
%  fn_sliderenhance(sliderhandle)

%% Description
%  allows a slider uicontrol to evaluate its callback during scrolling
%  (instead of only at the moment that the mouse if released)
%  just call 'fn_sliderenhance(sliderhandle)' once
% 
%  notes: 
%  - the function sets a 'WindowButtonDownFcn' property in the figure, hence
%    it will fail if this property is already set
%  - it also sets a 'DeleteFcn' to the slider and will fail if this property
%    is already set
%  - the functions needs to be called again each time the uicontrol is
%    resized
% 
%  See also fn_slider

%% Source
% Thomas Deneux
%
% Copyright 2009-2012
%
