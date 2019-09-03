function varargout = fn_mouse(varargin)
% function poly = fn_mouse([axes handle],'point|cross|poly|free|ellipse'[,msg])
% function [x y] = fn_mouse([axes handle],'point|cross'[,msg])
% function rect = fn_mouse([axes handle],'rect'[,msg])
% function [center axis e] = fn_mouse([axes handle],'ellipse'[,msg])
%---
% multi-functions function using mouse events
% mode defines action:
% 'point'       [default] get coordinates on mouse click
% 'cross'       get coordinates on mouse click - use cross pointer
% 'rect'        get a rectangle selection (format [xstart ystart xsize ysize])
% 'rectangle'   get a rectangle selection (format [x1 x2 x3 x4; y1 y2 y3 y4])
% 'poly'        polygone selection
% 'free'
% 'ellipse'
% options: (ex: 'rect+', 'poly-@:.25:')
% +     selection is drawn (all modes)
% -     use as first point the current point in axes (rect, poly, free, ellipse)
% @     plots open line instead of closed polygon (poly, free, ellipse)
% :num: interpolates line with one point every 'num' pixel (poly, free, ellipse)
%       for 'ellipse' mode, if :num; is not specified, output is a cell
%       array {center axis e} event in the case of only one outpout argument
% 
% See also fn_maskselect

% Thomas Deneux
% Copyright 2005-2012

% Input
i=1;
ha=[];
mode='';
msg = '';
while i<=nargin
    arg = varargin{i};
    if ishandle(arg), ha=arg;
    elseif ischar(arg)
        if isempty(mode), mode=arg; else msg=arg; end
    else error('bad argument')
    end
    i=i+1;
end
if isempty(mode), mode='point'; end
if isempty(ha), ha=gca; end
hf = get(ha,'parent');
while ~strcmp(get(hf,'type'),'figure'), hf = get(hf,'parent'); end
figure(hf)

% Extract parameters from mode definition
type = regexp(mode,'^(\w)+\>','match'); type = type{1};
buttonalreadypressed = ~isempty(strfind(mode,'-'));
showselection = strfind(mode,'+');
openline = strfind(mode,'@');
dointerp = strfind(mode,':');

switch type
    case 'point'
        SuspendCallbacks(ha)
        waitforbuttonpressmsg(ha,msg)
        RestoreCallbacks(ha)
        point = get(ha,'CurrentPoint');    % button down detected
        if nargout<=1, varargout={point(1,1:2)'};
        elseif nargout==2, varargout=num2cell(point(1,1:2));
        elseif nargout==3, varargout=num2cell(point(1,1:3));
        end
        if showselection
            oldnextplot=get(ha,'NextPlot'); set(ha,'NextPlot','add')
            plot(point(1,1),point(1,2),'+','parent',ha),
            set(ha,'NextPlot',oldnextplot)
        end
    case 'cross'
        SuspendCallbacks(ha)
        point = ginput(1);
        RestoreCallbacks(ha)
        if nargout<=1, varargout={point(1,:)'};
        elseif nargout==2, varargout=num2cell(point(1,1:2));
        elseif nargout==3, varargout=num2cell(point(1,1:3));
        end
        if showselection,
            oldnextplot=get(ha,'NextPlot'); set(ha,'NextPlot','add')
            plot(point(1,1),point(1,2),'+','parent',ha),
            set(ha,'NextPlot',oldnextplot)
        end
    case 'line'
        if ~buttonalreadypressed, waitforbuttonpressmsg(ha,msg), end
        p1 = get(ha,'CurrentPoint'); p1 = p1(1,1:2);
        hl(1) = line('xdata',p1([1 1]),'ydata',p1([2 2]),'parent',ha,'color','k','erasemode','xor');
        hl(2) = line('xdata',p1([1 1]),'ydata',p1([2 2]),'parent',ha,'color','b','erasemode','xor','linestyle','--');
        data = fn_buttonmotion({@drawline,ha,hl,p1},hf);
        delete(hl(2))
        if showselection
            set(hl(1),'color','y','erasemode','normal')
        else
            delete(hl(1))
        end
        if nargout<=1, varargout = {data};
        else varargout = {data(:,1) data(:,2)};
        end
    case {'rect','rectangle'}
        % if button has already been pressed, no more button will be
        % pressed, so it is not necessary to suspend callbacks
        SuspendCallbacks(ha)
        if ~buttonalreadypressed, waitforbuttonpressmsg(ha,msg), end
        p0 = get(ha,'currentpoint'); p0 = p0(1,1:2);
        hl(1) = line(p0(1),p0(2),'color','k','linestyle','-','parent',ha);
        hl(2) = line(p0(1),p0(2),'color','w','linestyle',':','parent',ha);
        rect = fn_buttonmotion({@drawrectangle,ha,hl,p0},hf);
        delete(hl)
        RestoreCallbacks(ha)
        if showselection,
            line(rect(1,[1:4 1]),rect(2,[1:4 1]),'color','k','parent',ha)
            line(rect(1,[1:4 1]),rect(2,[1:4 1]),'color','w','linestyle',':','parent',ha),
        end
        if nargout<=1
            if strcmp(type,'rectangle')
                varargout={rect};
            else % type is 'rect'
                cornera = [min(rect(1,:)); min(rect(2,:))];
                cornerb = [max(rect(1,:)); max(rect(2,:))];
                rect = [cornera cornerb-cornera];
                varargout = {rect};
            end
        end
    case 'poly'
        SuspendCallbacks(ha)
        if ~buttonalreadypressed, waitforbuttonpressmsg(ha,msg), end
        [xi,yi] = fn_getline(ha);                   % return figure units
        RestoreCallbacks(ha)
        x = [xi yi];
        if showselection,
            if openline, back=[]; else back=1; end
            oldnextplot=get(ha,'NextPlot'); set(ha,'NextPlot','add')
            plot(x([1:end back],1),x([1:end back],2),'k-'),
            plot(x([1:end back],1),x([1:end back],2),'w:'),
            set(ha,'NextPlot',oldnextplot)
        end
        if dointerp
            f = find(mode==':'); f = [f length(mode)+1];
            ds = str2num(mode(f(1)+1:f(2)-1));
            if ~openline, x(end+1,:)=x(1,:); end
            ni = size(x,1);
            L = zeros(ni-1,1);
            for i=2:ni, L(i) = L(i-1)+norm(x(i,:)-x(i-1,:)); end
            if ~isempty(L), x = interp1(L,x,0:ds:L(end)); end
        end
        if nargout==1, varargout={x'}; end
    case 'free'
        SuspendCallbacks(ha)
        if ~buttonalreadypressed, waitforbuttonpressmsg(ha,msg), end
        p = get(ha,'currentpoint');
        hl(1) = line(p(1,1),p(1,2),'color','k','linestyle','-','parent',ha);
        hl(2) = line(p(1,1),p(1,2),'color','w','linestyle',':','parent',ha);
        fn_buttonmotion({@freeform,ha,hl},hf)
        RestoreCallbacks(ha)    
        x = [get(hl(1),'xdata')' get(hl(2),'ydata')'];
        delete(hl)
        if showselection,
            if openline, back=[]; else back=1; end
            oldnextplot=get(ha,'NextPlot'); set(ha,'NextPlot','add')
            plot(x([1:end back],1),x([1:end back],2),'k-','parent',ha),
            plot(x([1:end back],1),x([1:end back],2),'w:','parent',ha),
            set(ha,'NextPlot',oldnextplot)
        end
        if dointerp
            f = find(mode==':'); f = [f length(mode)+1];
            ds = str2num(mode(f(1)+1:f(2)-1));
            if ~openline, x(end+1,:)=x(1,:); end
            ni = size(x,1);
            L = zeros(ni-1,1);
            for i=2:ni, L(i) = L(i-1)+norm(x(i,:)-x(i-1,:)); end
            if ~isempty(L), x = interp1(L,x,0:ds:L(end)); end
        end
        if nargout==1, varargout={x'}; end
    case 'ellipse'
        SuspendCallbacks(ha)
        if ~buttonalreadypressed, waitforbuttonpressmsg(ha,msg), end
        p = get(ha,'currentpoint');
        hl(1) = line(p(1,1),p(1,2),'color','k','linestyle','-','parent',ha);
        hl(2) = line(p(1,1),p(1,2),'color','w','linestyle',':','parent',ha);
        pflag = pointer('init');
        % circle
        fn_buttonmotion({@drawellipse,ha,hl,pflag},hf)
        if strcmp(getvalue(pflag),'width')
            % change eccentricity -> ellipse
            fn_buttonmotion({@drawellipse,ha,hl,pflag},hf)
        end
        RestoreCallbacks(ha)   
        x = [get(hl(1),'xdata')' get(hl(1),'ydata')'];
        ax = getappdata(hl(1),'axis');
        axis = (ax(:,2)-ax(:,1))/2;
        center = getappdata(hl(1),'start') + axis;
        e = getappdata(hl(1),'excentricity');
        delete(hl)
        if showselection,
            if openline, back=[]; else back=1; end
            oldnextplot=get(ha,'NextPlot'); set(ha,'NextPlot','add')
            plot(x([1:end back],1),x([1:end back],2),'k-','parent',ha),
            plot(x([1:end back],1),x([1:end back],2),'w:','parent',ha),
            set(ha,'NextPlot',oldnextplot)
        end
        if dointerp
            f = find(mode==':'); f = [f length(mode)+1];
            ds = str2num(mode(f(1)+1:f(2)-1));
            if ~openline, x(end+1,:)=x(1,:); end
            ni = size(x,1);
            L = zeros(ni-1,1);
            for i=2:ni, L(i) = L(i-1)+norm(x(i,:)-x(i-1,:)); end
            if ~isempty(L), x = interp1(L,x,0:ds:L(end)); end
        end
        switch nargout
            case 1
                if dointerp
                    varargout = {x'};
                else
                    varargout = {{center axis e}};
                end
            case 3
                varargout = {center axis e};
        end
    otherwise
        error('unknown type ''%s''',type)
end

%-------------------------------------------------
function SuspendCallbacks(ha)
% se pr�munir des callbacks divers et vari�s

hf = get(ha,'parent'); while ~strcmp(get(hf,'type'),'figure'), hf = get(hf,'parent'); end
setappdata(ha,'uistate',guisuspend(hf))
setappdata(ha,'oldtag',get(ha,'Tag'))
set(ha,'Tag','fn_mouse') % pour bloquer fn_imvalue !

%-------------------------------------------------
function RestoreCallbacks(ha)
% r�tablissement des callbacks avant les affichages

set(ha,'Tag',getappdata(ha,'oldtag'))
rmappdata(ha,'oldtag')
guirestore(getappdata(ha,'uistate'))

%-------------------------------------------------
function state = guisuspend(hf)

state.hf        = hf;
state.obj       = findobj(hf);
state.hittest   = get(state.obj,'hittest');
state.buttonmotionfcn   = get(hf,'windowbuttonmotionfcn');
state.buttondownfcn     = get(hf,'windowbuttondownfcn');
state.buttonupfcn       = get(hf,'windowbuttonupfcn');
state.keydownfcn        = get(hf,'keypressfcn');
try state.keyupfcn = get(hf,'keyreleasefcn'); end

set(state.obj,'hittest','off')
set(hf,'hittest','on','windowbuttonmotionfcn','', ...
    'windowbuttondownfcn','','windowbuttonupfcn','', ...
    'keypressfcn','')
try set(hf,'keyreleasefcn',''), end

%-------------------------------------------------
function guirestore(state)

for k=1:length(state.obj)
    set(state.obj(k),'hittest',state.hittest{k});
end
hf = state.hf;
set(hf,'windowbuttonmotionfcn',state.buttonmotionfcn);
set(hf,'windowbuttondownfcn',state.buttondownfcn);
set(hf,'windowbuttonupfcn',state.buttonupfcn);
set(hf,'keypressfcn',state.keydownfcn);
try set(hf,'keyreleasefcn',state.keyupfcn); end

%-------------------------------------------------
function data=drawline(ha,hl,p1)

p2 = get(ha,'currentpoint');
data = [p1(:) p2(1,1:2)'];
set(hl,'xdata',data(1,:),'ydata',data(2,:))

%-------------------------------------------------
function freeform(ha,hl)

p = get(ha,'currentpoint');
xdata = get(hl(1),'xdata'); xdata(end+1) = p(1,1);
ydata = get(hl(1),'ydata'); ydata(end+1) = p(1,2);
set(hl,'xdata',xdata,'ydata',ydata)

%-------------------------------------------------
function rect = drawrectangle(ha,hl,p0)

p = get(ha,'currentpoint'); p = p(1,1:2);
xdata = [p([1 1]) p0([1 1]) p(1)];
ydata = [p(2) p0([2 2]) p([2 2])];
rect = [xdata(1:4); ydata(1:4)];
set(hl,'xdata',xdata,'ydata',ydata)

%-------------------------------------------------
function drawellipse(ha,hl,pflag)

p = get(ha,'currentpoint');
switch getvalue(pflag)
    case 'init'
        xdata = get(hl(1),'xdata');
        ydata = get(hl(1),'ydata');
        setappdata(hl(1),'start',[xdata; ydata])
        setvalue(pflag,'axis')
        set(get(ha,'parent'),'windowbuttondownfcn',@(hf,evnt)setvalue(pflag,'click'))
        drawellipse(ha,hl,pflag)
        return
    case 'axis'
        ax = [getappdata(hl(1),'start') p(1,1:2)'];
        u = (ax(:,2)-ax(:,1))/2;
        e = .999; % make the 2 eigenvalues different so that the 2 eigenvectors will not be changed once automatically recalculated
    case 'click'
        ax = getappdata(hl(1),'axis');
        u = (ax(:,2)-ax(:,1))/2;
        v = [u(2); -u(1)];
        p = ax(:,1) + u + v;
        p0 = fn_coordinates(ha,'a2s',p,'position');
        set(0,'pointerlocation',p0);
        setvalue(pflag,'width')
        return
    case 'width'
        ax = getappdata(hl(1),'axis');
        u = (ax(:,2)-ax(:,1))/2;
        v = [u(2); -u(1)];
        normu2 = sum(u.^2);
        o = mean(ax,2);
        op = p(1,1:2)'-o;
        uc = sum(op.*u)/normu2;
        vc = sum(op.*v)/normu2;
        e = abs(vc / (sin(acos(uc))));
end
v = [u(2) -u(1)];
t = 0:.05:1;
udata = (1-cos(2*pi*t));
vdata = e*sin(2*pi*t);
xdata = ax(1,1) + u(1)*udata + v(1)*vdata;
ydata = ax(2,1) + u(2)*udata + v(2)*vdata;
set(hl,'xdata',xdata,'ydata',ydata)
setappdata(hl(1),'axis',ax)
setappdata(hl(1),'excentricity',e)

%---
function waitforbuttonpressmsg(ha,msg)

hf = get(ha,'parent');

%if isempty(msg), waitfor(hf,'windowbuttondownfcn',''), return, end

p = get(ha,'currentpoint'); p=p(1,1:2);
dd = fn_coordinates(ha,'b2a',[9 -12],'vector')';
t = text('parent',ha,'string',msg, ...
    'fontsize',8,'erasemode','xor', ...
    'position',p(1:2)+dd);
set(hf,'windowbuttonmotionfcn',@(f,evnt)movesub(ha,t,dd), ...
    'windowbuttondownfcn', ...
    @(f,evnt)set(hf,'windowbuttonmotionfcn','','windowbuttondownfcn',''))
waitfor(hf,'windowbuttondownfcn','')
delete(t)

%---
function movesub(ha,t,dd)

p = get(ha,'currentpoint'); p=p(1,1:2);
set(t,'position',p(1:2)+dd)

