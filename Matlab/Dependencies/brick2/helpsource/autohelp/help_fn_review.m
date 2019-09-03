%% fn_review

%% Syntax
%  fn_review(x1,x2,...[,command])
%  fn_review(x[,command])

%% Description
%  opens a new figure and displays one of the datas xk (change k by pressing
%  arrows); if there is only one data x, switches between subdata according
%  to the following rules:
%  - if x is a structure or an object, uses xk = x(k) [k can be multidimensional index]
%  - if x is a cell array, uses xk = x{k} [idem]
%  - if x is an array, uses xk = x(:,..,k), operating on the last dimension
% 
%  if there is a command argument, executes the custom command instead
%  of the default
%  [technical note: the custom command can be either
%   - a character array to be evaluated in base workspace using variable 'x'
%   - or a function handle - fuction should have 1 argument]
%  
%  the default command is @showrew, file fn_review_showres.m can be edited to change
%  the behaviour
% 
%  see also fn_review_showres

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
