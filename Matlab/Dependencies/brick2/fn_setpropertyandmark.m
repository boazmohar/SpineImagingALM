function fn_setpropertyandmark(C,propertyname,hu,value) %#ok<INUSL>
% function fn_setpropertyandmark(C,propertyname,hu,value)
%---
% This function is useful for setting the callback of a boolean uibutton,
% or of a menu with a checkmark which reflects the boolean value of a
% specific property of a given object.
% Once invoked, it set both the property value of the object, and the
% control display.
% 
% Input:
% - C               object
% - propertyname    property name; can be of the form 'prop.subfield1.etc'
% - hu              control or uimenu handle
% - value           boolean, or 'switch' to inverse the value

% Thomas Deneux
% Copyright 2007-2012

if ischar(value) 
    if ~strcmp(value,'switch'), error argument, end
    b = ~eval(['C.' propertyname]);
else
    b = value;
end

% object property
valstr = fn_switch(b,'true','false');
eval(['C.' propertyname '=' valstr])
 
% control display
switch get(hu,'type')
    case 'uimenu'
        set(hu,'checked',fn_switch(b,'on','off'))
    case 'uicontrol'
        set(hu,'value',b)
    otherwise
        error('wrong handle')
end
