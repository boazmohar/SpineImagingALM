function filename = fn_savefile(varargin)
% function filename = fn_savefile(varargin)
%--
% synonyme de "filename = fn_getfile('SAVE',varargin)"
% 
% See also fn_getfile

% Thomas Deneux
% Copyright 2003-2012

filename = fn_getfile('SAVE',varargin{:});