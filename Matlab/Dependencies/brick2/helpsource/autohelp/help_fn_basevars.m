%% fn_basevars

%% Syntax
%  fn_basevars(['get|send',]name1,name2)

%% Description
%  load base workspace in caller workspace or vice-versa (if 'send' flag)
%  usefull when debugging a function and trying to use it as a script
%  load variables whose name are given as arguments, or the whole base 
%  workspace with noargument 
%  Performs inverse operation if first argument is numeric

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
