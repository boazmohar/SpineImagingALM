function out = fn_fileparts(f,flag)
% function out = fn_filename(fname,flag)
%---
%
% Input:
% - fname       file name
% - flag        one of the possible flags below
% 
% 'path'    returns the path (fn_filename('a/b/','path') returns 'a')
% 'base'    removes the path and the extension (fn_filename('a/b.c','base')
%           returns 'b') 
% 'ext'     returns the extension (fn_filename('a/b.c','ext') returns '.c')
% 'name'    removes the path (fn_filename('a/b/','name') returns 'b',
%           fn_filename('a/b.c','name') returns 'b.c') 
% 'noext'   removes the extention (fn_filename('a/b.c','noext') returns
%           'a/b') 
% ''        only removes trailing '/' characters

% Thomas Deneux
% Copyright 2003-2012

while f(end)=='/', f(end)=[]; end
[p b e] = fileparts(f);
switch flag
    case ''
        out = f;
    case 'path'
        out = p;
    case 'base'
        out = b;
    case 'ext'
        out = e;
    case 'name'
        out = fullfile(p,b);
    case 'noext'
        out = fullfile(p,b);
    otherwise
        error('unknown flag ''%s''',flag)
end
