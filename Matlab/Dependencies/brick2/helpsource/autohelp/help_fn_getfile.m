%% fn_getfile

%% Syntax
%  filename = fn_getfile(['READ|SAVE',][filter[,title]])
%  filename = fn_getfile('DIR',title)
%  rep = fn_getfile('REP')
%  fn_getfile('REP',rep)

%% Description
%  returns file name with full path
%  fn_getfile('REP') returns the current directory
%  fn_getfile('REP',rep) sets the current directory
%  as a utility, fn_getfile('cd') is going to the current directory 
%  
%  See also fn_savefile

%% Source
% Thomas Deneux
%
% Copyright 2003-2012
%
