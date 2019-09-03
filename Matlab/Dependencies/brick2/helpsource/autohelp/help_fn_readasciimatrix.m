%% fn_readasciimatrix

%% Syntax
%  a=fn_readasciimatrix(filename,nheaders)

%% Description
%  
%  load a matrix from an ascii file with format like :
%    #blah blah
%    et reblah 45 blah
%    1.2 3 5.5
%    4   0 0 
%    3   3 14
%  empty matrix is returned if file does not exist or if there is no numeric
%  content at any line begin
%  there is no verification that each row is the same length
% 
%  if 'nheader' is specified, first nheaders lines are skipped anyway
% 
%  See also fn_saveasciimatrix, fn_readdatlabview, fn_readtext, fn_readbin

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
