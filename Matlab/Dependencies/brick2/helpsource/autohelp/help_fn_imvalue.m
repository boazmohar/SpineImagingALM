%% fn_imvalue

%% Description
%  fn_imvalue [image] [xy|xonly]
%  fn_imvalue clean
%  fn_imvalue end
%  fn_imvalue('chgx'|'chgy'|'chgxy',newax[,ha])
%  fn_imvalue('register'[,ha][,command])
%  fn_imvalue('unregister'[,ha])
%  fn_imvalue demo
%  Link images with same dimensions via a crosshair pointer, prints
%  values, and enable common zooming and translations.
%  Link plots with same x-axis and enable common x-zooming and
%  x-translations.
%  To prevent conflicts with other programs, it does nothing in windows
%  that contain objects with a 'Tag' or 'ButtonDown' property already set.
% 
%  When an axes is clicked, its 'UserData' property is changed, such that
%  user can add additional callbacks to the same axes by using a listener
%  ('addlistener(ha,'UserData','PostSet',fun_handle)') - don't use the axes
%  'ButtonDownFcn' property as it is already used by fn_imvalue.

%% Source
% Thomas Deneux
%
% Copyright 2003-2012
%
