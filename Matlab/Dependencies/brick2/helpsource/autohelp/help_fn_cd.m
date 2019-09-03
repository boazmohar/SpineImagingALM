%% fn_cd

%% Syntax
%  rep = fn_cd(flag,relpath)
%  fn_cd('edit')
% This  allows user to rapidly change the current directory or get

%% Description
%  the full path to some files. 
% 
%  - Type 'fn_cd edit' to launch a GUI interface that enables you to define
%    your own absolute and relative paths (for example the full path to your
%    home, to which will be associated the flag 'home')
%  - Then type 'fn_cd home' to set Matlab current directory to your home.
%  - Type 'a = fn_cd('home');' to get the full path to your home (note that
%    this will not change Matlab current directory
%  - Additional arguments to fn_cd can be used to access subdirectories or
%    files inside a flagged folder. For example, type 'fn_cd home dir1 dir2'
%    to set Matlab current directory to the subdirectory dir1/dir2 inside
%    your home directory, or type 'a = fn_cd('home','dir1/myfile');' to get
%    the full path to the file 'myfile' inside the subdirectory 'dir1'
%    inside your home directory.

%% Source
% Thomas Deneux
%
% Copyright 2002-2012
%
