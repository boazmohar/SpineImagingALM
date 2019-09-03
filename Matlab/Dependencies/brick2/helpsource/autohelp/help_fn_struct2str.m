%% fn_struct2str

%% Syntax
%  str = fn_struct2str(par[,parname[,'cell'])

%% Description
%  Input:
%  - par             structure
%  - parname         name (use inputname(1) or 'x' if not specified)
%  - 'cell' flag     returns a cell array of srings instead of a character
%                    array
% 
%  Output:
%  - str         character array which fills the structure with values
% 
%  ex: par = struct('a',2,'b','hello'); fn_struct2str(par,'x')
%    returns: 
%      x.a = 2;
%      x.b = hello
% 
%  See also fn_str2struct, fn_structedit

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
