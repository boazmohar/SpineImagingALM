function fn_evalcallback(fun,hobj,evnt)
% function fn_evalcallback(fun,hobj,evnt)

% Thomas Deneux
% Copyright 2007-2012

if nargin<2, hobj = []; end
if nargin<3, evnt = []; end

switch class(fun)
    case 'char'
        evalin('base',fun)
    case 'function_handle'
        feval(fun,hobj,evnt)
    case 'cell'
        feval(fun{1},hobj,evnt,fun{2:end})
    otherwise
        error('graphic object callback cannot be of class ''%s''',class(fun))
end
            