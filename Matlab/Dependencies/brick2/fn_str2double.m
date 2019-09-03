function x = fn_str2double(x)
% function x = fn_str2double(x)
%---
% scan character array x to read a double... only if is really a character
% array!!
%
% See also fn_num2str

% Thomas Deneux
% Copyright 2007-2012

if ischar(x)
    x = str2double(x);
end