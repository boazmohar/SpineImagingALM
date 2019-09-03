%% fn_supercontrol

%% Syntax
%  X = fn_supercontrol([hp,]specs[,callback[,x0]])
%  specs = fn_supercontrol.examplespecs

%% Description
%  Input:
%  - hp          uipanel handle
%  - specs       structure with fields:
%                . name      string
%                . controls  structure with fields (* is mandatory):
%                            style*  popupmenu, edit, checkbox,
%                                    pushbutton or stepper
%                            string  the 'string' property of the
%                                    control
%                            length* relative width occupied by the
%                                    control + its label if any
%                            default* default value (type depends on
%                                    control style)
%                            label   a label placed on the left
%                            labellength     relative width occupied by
%                                    the label (set to 0 if no label)
%                            callback        for push button only:
%                                    function with prototype 
%                                    valuek = fun(value) where value is
%                                    a cell array and valuek an element
%                                    of this array [MORE DOC NEEDED]
%                            more    more properties stored in a cell
%                                    array with successive pairs of
%                                    property names/values
%  - callback    function to be executed when control values are
%                changed by user, with prototype @(x)fun(x),
%                where x is X.x (see below)
%  - x0          initial value (see below comments on X.x)
% 
%  Output:
%  - X           fn_supercontrol object; X.x is a structure that
%                stores the values, with fields:
%                .name       string
%                .active     logical
%                .value      cell array with values (one per
%                            control in the specs of the same
%                            name)
% 
%  One can change the values either by user action (acionning
%  the controls) or by setting the value of X.x.
% 
%  Special notes:
%  - chgactive flag: The callback might be invoked but the
%  "active" part of the value has not changed (for example, when
%  a new, inactive line has been creadted). The property
%  X.activechg says whether this active part has changed or not.
%  - 'edit': if a given control has the style 'edit', but its
%  data is set to numeric value instead of a string, the data is
%  stored in the 'userdata' property of the control, and the
%  string 'userdata' is diplayed.
%  - 'pushbutton': when pressing the button, the callback is executed
%  and changes the value accordingly
% 
%  See also fn_control

%% Source
% Thomas Deneux
%
% Copyright 2010-2012
%
