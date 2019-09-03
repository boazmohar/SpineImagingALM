%% fn_setpropertyandmark

%% Syntax
%  fn_setpropertyandmark(C,propertyname,hu,value)
% This  is useful for setting the callback of a boolean uibutton,

%% Description
%  or of a menu with a checkmark which reflects the boolean value of a
%  specific property of a given object.
%  Once invoked, it set both the property value of the object, and the
%  control display.
%  
%  Input:
%  - C               object
%  - propertyname    property name; can be of the form 'prop.subfield1.etc'
%  - hu              control or uimenu handle
%  - value           boolean, or 'switch' to inverse the value

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
