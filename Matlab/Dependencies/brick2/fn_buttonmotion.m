function varargout = fn_buttonmotion(fun,hf,varargin)
% function varargout = fn_buttonmotion(fun[,hf])
% function fn_buttonmotion('demo')
%---
% utility for executing a task while mouse button is pressed, and execute
% it a last time when button is realeased
% 
% TODO: add custom queuing system
%
% See also fn_moveobject

% Thomas Deneux
% Copyright 2007-2012

% Old version
if nargin>2 || (nargin>1 && ~ishandle(hf))
    error('Syntax for fn_buttonmotion function has changed, please replace fn_buttonmotion(@myfun,arg1,arg2,...) by fn_buttonmotion({@myfun,arg1,arg2,..})')
end

disp_if_debug('entering fn_buttonmotion function')

% Input
if ischar(fun) && strcmp(fun,'demo')
    demo()
    return
end
if nargin<2
    hf = gcf;
end

% Already an existing callback!
if ~isempty(get(hf,'WindowButtonMotionFcn')) ...
        ||  ~isempty(get(hf,'WindowButtonUpFcn'))
    disp_if_debug('Problem! WindowButtonMotionFcn already set to: ',get(hf,'WindowButtonMotionFcn'))
    disp_if_debug('Problem! WindowButtonUpFcn already set to:     ',get(hf,'WindowButtonUpFcn'))
    terminate(hf)
    warning('a WindowButtonMotionFcn, DownFcn or UpFcn property is already set to figure') %#ok<WNTAG>
end
okdown = isempty(get(hf,'WindowButtonDownFcn'));

% Motion
disp_if_debug('current WindowButtonMotionFcn is: ',get(hf,'WindowButtonMotionFcn'),', setting now to @motionexec')
set(hf,'WindowButtonMotionFcn',{@motionexec 'move' fun}, ...
    'WindowButtonUpFcn',{@motionexec 'stop'})
if okdown, set(hf,'WindowButtonDownFcn',{@motionexec 'stop'}), end
setappdata(hf,'fn_buttonmotion_scrolling',true)
disp_if_debug('waiting for motion end')
waitfor(hf,'WindowButtonMotionFcn','')
if okdown, set(hf,'WindowButtonDownFcn',''), end

% Button release
if ~ishandle(hf), return, end % figure has been closed in the mean while
setappdata(hf,'fn_buttonmotion_scrolling',false)
[varargout{1:nargout}] = exec(fun);
rmappdata(hf,'fn_buttonmotion_scrolling') 

%---
function motionexec(hf,evnt,actionflag,fun)

% start execution
debugstr = [actionflag ' ' num2str(floor(100*rand))];
disp_if_debug(['start ' debugstr])

% custom queuing/canceling system
if getappdata(hf,'fn_buttonmotion_busy')
    switch actionflag
        case 'move' % cancel
            disp_if_debug(['rejct ' debugstr])
            return
        case 'stop' % queue
            setappdata(hf,'fn_buttonmotion_queue',debugstr)
            disp_if_debug(['queue ' debugstr])
            return
    end
end
setappdata(hf,'fn_buttonmotion_busy',true)

% stop (not queued)
if strcmp(actionflag,'stop')
    set(hf,'WindowButtonMotionFcn','', ...
        'WindowButtonUpFcn','')
    rmappdata(hf,'fn_buttonmotion_busy')
    disp_if_debug(['end   ' debugstr])
    return
end

% evaluate function
disp_if_debug(['exec  ' debugstr])
try exec(fun); catch ME, disp hello, terminate(hf), rethrow(ME), end

% end of current execution
disp_if_debug(['end   ' debugstr])
setappdata(hf,'fn_buttonmotion_busy',false)

% stop (queued)
debugstr = getappdata(hf,'fn_buttonmotion_queue');
if ~isempty(debugstr)
    set(hf,'WindowButtonMotionFcn','', ...
        'WindowButtonUpFcn','')
    rmappdata(hf,'fn_buttonmotion_busy')
    rmappdata(hf,'fn_buttonmotion_queue')
    disp_if_debug(['exec  ' debugstr])
end

%---
function varargout = exec(fun)

if ischar(fun)
    evalin('base',fun);
elseif isa(fun,'function_handle')
    [varargout{1:nargout}] = feval(fun);
elseif iscell(fun)
    [varargout{1:nargout}] = feval(fun{:});
else
    error bad
end

%---
function terminate(hf)

set(hf,'WindowButtonMotionFcn','', ...
    'WindowButtonUpFcn','')
try rmappdata(hf,'fn_buttonmotion_busy'), end %#ok<*TRYNC>
try rmappdata(hf,'fn_buttonmotion_queue'), end

%---
function disp_if_debug(varargin)

% str = [];
% for k=1:nargin
%     x = varargin{k};
%     if iscell(x)
%         if isa(x{1},'function_handle')
%             x = func2str(x{1});
%         else
%             error('don''t know how to display cell array')
%         end
%     end
%     str = [str x]; %#ok<AGROW>
% end
% disp(str)

%---
function demo

C = {'figure(1), clf'
    'ht = uicontrol(''style'',''text'');'
    'fun = ''p=get(1,''''currentpoint''''); p=p(1,1:2); set(ht,''''pos'''',[p 60 20],''''string'''',num2str(p))'';'
    'set(1,''buttondownfcn'',@(u,evnt)fn_buttonmotion(fun))'};
fn_dispandexec(C)


