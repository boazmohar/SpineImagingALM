function x = fn_input(name,varargin)
% function x = fn_input(name[,defaultval][,min,max])
%---
% this is a small wrapper for function fn_structedit, to prompt user for
% a single value; normally if defaultval, min and max are integers, x will
% also be integer
% 
% See also fn_structedit, fn_control, fn_reallydlg

% Thomas Deneux
% Copyright 2007-2012

% input
if nargin<1, name = 'x'; end
switch nargin
    case {0,1}
        defaultval = 0; 
        min = 0; 
        max = 5;
    case 2
        defaultval = varargin{1};
        if defaultval<0
            min = 2*defaultval;
            max = 2*abs(defaultval);
        elseif defaultval<100
            min = 0;
            max = 2*defaultval;
        else
            min = defaultval/10;
            max = defaultval*10;
        end
    case 3
        [min max] = deal(varargin{:});
        defaultval = min;
    case 4
        [defaultval min max] = deal(varargin{:});
end

% specifications
if length(defaultval)>1 || nargin<3
    spec = class(defaultval);
elseif min>0 &&  max/min>100
    % logarithmic scale
    min = log10(min);
    max = log10(max);
    spec = ['logslider ' num2str(min) ' ' num2str(max) ' .1 %.2g'];
elseif ~mod(defaultval,1) && ~mod(min,1) && ~mod(max,1)
    spec = ['stepper 1 ' num2str(min) ' ' num2str(max)];
else
    % slider, min, max
    spec = ['slider ' num2str(min) ' ' num2str(max)];
end

% call to fn_structedit
s = struct(name,defaultval);
spec = struct(name,spec);
s = fn_structedit(s,spec);

% output
if isempty(s), x=[]; else x=s.(name); end
            
