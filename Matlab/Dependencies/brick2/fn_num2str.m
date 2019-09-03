function s = fn_num2str(s,n)
% function s = fn_num2str(s)
%---
% convert numerical value s into a string representation... unless s is
% already a character array!
%
% See also fn_str2double

% Thomas Deneux
% Copyright 2007-2012

% overloaded function
if nargin==2
    disp('fn_num2str_old is desuete')
    s = fn_num2str_old(s,n);
    return
end

% this function
if ~ischar(s)
    s = num2str(s);
end