%% fn_switch

%% Syntax
%  y = fn_switch(x,y_true,[x2,y_true2,[...]]y_false)
%  y = fn_switch(x,case1,y1,case2,y2,..,casen,yn[,ydefault])
%  y = fn_switch(true|false)
%  y = fn_switch('on|off'[',switch'])
% the 2 first cases are general prototypes: the s recognize which

%% Description
%  to use according to whether x is scalar and logical
%  MAKE SURE THAT X IS SCALAR AND LOGICAL IF YOU WANT TO USE THE FIRST FORM!
%  the 2 other cases are specialized shortcuts to convert logical values
%  true or false in the string 'on' or 'off'
%  'case1', 'case2', etc.. can be any Matlab variable to which x is compared
%  but if they are cell arrays, x is compared to each of their elements

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
