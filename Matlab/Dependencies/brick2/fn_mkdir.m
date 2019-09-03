function y = fn_mkdir(D)
% function y = fn_mkdir(D)
% creates directory D if it does not exist
% D can be an asolute or relative path
% parent directories are also created if necessary

% Thomas Deneux
% Copyright 2004-2012

spwd=pwd;

Dpath = fn_fileparts(D,'path');
if isempty(Dpath)
    Dpath = pwd;
end
Dend = fn_fileparts(D,'name');

try 
    
    cd(Dpath)
    warning off MATLAB:MKDIR:DirectoryExists
    mkdir(Dend)
    warning on MATLAB:MKDIR:DirectoryExists
    
catch
    
    fn_mkdir(Dpath)
    fn_mkdir(D)
    
end

cd(spwd)