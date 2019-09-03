function b = fn_reallydlg(varargin)
% function b = fn_reallydlg(line1,line2,..)
%--
% returns true if 'Yes' has been answered
% argument can be one or several string or cell array of strings
% 
% See also fn_input, fn_structedit

% Thomas Deneux
% Copyright 2009-2012

for i=1:nargin
    if ~iscell(varargin{i}), varargin{i}={varargin{i}}; end
end
question = [varargin{:}];
answer = questdlg(question,'warning','Yes','No','No');
b = strcmp(answer,'Yes');
