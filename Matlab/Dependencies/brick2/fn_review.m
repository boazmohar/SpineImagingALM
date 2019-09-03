function fn_review(varargin)
% function fn_review(x1,x2,...[,command])
% function fn_review(x[,command])
%---
% opens a new figure and displays one of the datas xk (change k by pressing
% arrows); if there is only one data x, switches between subdata according
% to the following rules:
% - if x is a structure or an object, uses xk = x(k) [k can be multidimensional index]
% - if x is a cell array, uses xk = x{k} [idem]
% - if x is an array, uses xk = x(:,..,k), operating on the last dimension
%
% if there is a command argument, executes the custom command instead
% of the default
% [technical note: the custom command can be either
%  - a character array to be evaluated in base workspace using variable 'x'
%  - or a function handle - fuction should have 1 argument]
% 
% the default command is @showrew, file fn_review_showres.m can be edited to change
% the behaviour
%
% see also fn_review_showres

% Thomas Deneux
% Copyright 2005-2012

if nargin==0, help fn_review, return, end

% Input
command = varargin{end};
if ischar(command) || isa(command,'function_handle')
    narg = nargin-1;
else
    command = 'default';
    narg = nargin;
end
if narg>1
    X = varargin(1:narg);
else
    X = varargin{1};
    if iscell(X)
        % X is already a good boy
    elseif isstruct(X) || isobject(X)
        X = num2obj(X);
    elseif isnumeric(X)
        X = num2cell(X,1:ndims(X)-1);
    else
        error('single data argument must be a cell array or struct array')
    end
end
    
[ni nj nk] = size(X); % important to do like this in case ndims(X)>3
s = [ni nj nk];
singleelem = false;
switch sum((s>1).*[1 2 4])
    case 0
        disp('fn_review: only one element')
        singleelem = true;
        d = find(s>1,1,'first');
        dims = [1 1 1];
        show = 1;
        fact = [1 1 1];
    case {1,2,4}    % vector
        d = find(s>1,1,'first');
        dims = [d d d];
        show = d;
        fact = [10 1 100];
    case 7          % 3D
        dims = [1 2 3];
        show = [1 2 3];
        fact = [1 1 1];
    case 3          % 2D
        if s(1)>=s(2)
            dims = [1 2 1];
        else
            dims = [1 2 2];
        end
        show = [1 2];
        fact = [1 1 10];
    otherwise
        error('stupid size')
end

% init
hf = figure(797);
set(hf,'numbertitle','off','name','fn_review')
clf(hf), fn_figmenu
info = struct('X',{X},'idx',[1 1 1], ...
    's',s,'dims',dims,'show',show,'fact',fact,'command',command, ...
    'ha',axes('parent',hf));
setappdata(hf,'fn_review',info)

if ~singleelem
    set(hf,'WindowKeyPressFcn',@keypress,'WindowScrollWheelFcn',@keypress)
end

evalcommand(info)

%---
function keypress(hf,evnt)

info = getappdata(hf,'fn_review');
if isfield(evnt,'VerticalScrollCount')
    sgn = evnt.VerticalScrollCount;
    fact = 1;
    d = 1;
else
    switch evnt.Key
        case {'leftarrow','space','rightarrow'}
            sgn = 1-2*strcmp(evnt.Key,'leftarrow');
            d = 2;
        case {'uparrow','downarrow'}
            sgn = 1-2*strcmp(evnt.Key,'uparrow');
            d = 1;
        case {'pageup','pagedown'}
            sgn = 1-2*strcmp(evnt.Key,'pageup');
            d = 3;
        case 'c'
            if strcmp(evnt.Modifier,'control'), close(hf), end
            return
        otherwise
            return
    end
    fact = info.fact(d);
end
d = info.dims(d);
info.idx(d) = 1+mod((info.idx(d)-1)+sgn*fact,info.s(d));

setappdata(hf,'fn_review',info)

evalcommand(info);

%---
function evalcommand(info)

try
    idx = info.idx;
    x = info.X{idx(1),idx(2),idx(3)};
    if isa(info.command,'function_handle')
        info.command(x);
    elseif strcmp(info.command,'default')
        fn_review_showres(x,info)
    else
        assignin('base','x',x)
        evalin('base',info.command)
        axis(info.ha,'tight')
    end
    title(info.ha,num2str(info.idx(info.show)))
catch ME
    disp(['Error in fn_review display: ' ME.message])
    title(info.ha,[num2str(info.idx(info.show)) ' [ERROR OCCURED]'])
end
