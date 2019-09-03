function x = fn_round(x,y)
% function x = fn_round(x[,y])
%---
% eliminate floating point error
% rounds to the nearest multiple of y if a second argument is given

% Thomas Deneux
% Copyright 2009-2012

if nargin==0, help fn_round, end

if nargin==2
    x = round(x/y)*y;
end

switch class(x)
    case 'double'
        x = str2num(num2str(x,'%.14g ')); %#ok<ST2NM>
    case 'single'
        x = str2num(num2str(x,'%.7g ')); %#ok<ST2NM>
    otherwise
        disp('strange usage of fn_round')
end


