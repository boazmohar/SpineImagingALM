function [pos cmd] = fn_framedesign(grob,pos,resetflag)
% function [pos cmd] = fn_framedesign([grob[,pos|cmd[,resetflag]]])
%---
% utility for positioning frames in a figure
%
% Input:
% - grob        structure with fields the names of graphic objects
%               considered and values their handles; their must be one
%               field 'hf' for the main figure
% - pos         structure with fields the names of graphic objects and
%               values their positions; fields are allowed to differ from
%               those of 'grob'
% - cmd         string defining pos (should be something like
%               'pos.hf=[..]; ...')
% - resetflag   true [default] or false - user resize?
%
% Output:
% - pos         as above
% - cmd         as above

% Thomas Deneux
% Copyright 2007-2012

% Input - make grob a structure with a 'hf' field
if nargin<1, grob = gcf; end
if nargin<2, pos = []; end
if nargin<3, resetflag = true; end
if length(grob)==1 && ishandle(grob) && strcmp(get(grob,'type'),'figure')
    hf = grob;
    obj = get(hf,'children')';
    obj = obj(isprop(obj,'units'));
    grob = struct('hf',hf,'obj',obj);
elseif all(ishandle(grob))
    obj = struct2array(grob);
    i = find(strcmp(get(obj,'type'),'figure'),1,'first');
    if isempty(i), obj(end+1)=get(obj(1),'parent'); i=length(obj); end
    hf = obj(i); obj = obj(setdiff(1:length(obj),i));
    grob = struct('hf',hf,'obj',obj); 
elseif ~isfield(grob,'hf')
    F=fieldnames(grob); grob.hf = get(grob.(F{1})(1),'parent'); 
end
figure(grob.hf)

% Set positions according to 'pos'
newgrob = true;
if nargin>1
    if ischar(pos), eval(pos), end
    newgrob = false;
    F = fieldnames(grob);
    for k=1:length(F)
        f = F{k};
        if ~isfield(pos,f) || size(pos.(f),1)~=numel(grob.(f))
            newgrob = true;
        end
        if isfield(pos,f)
            for i=1:min(size(pos.(f),1),numel(grob.(f)))
                set(grob.(f)(i),'position',pos.(f)(i,:))
            end
        end
    end
end

% Correct figure position however to make it visible if necessary
screensize = get(0,'screenSize');
if ~isempty(pos) && isfield(pos,'hf') && ...
        (pos.hf(1)+pos.hf(3)>screensize(3) || pos.hf(2)+pos.hf(4)>screensize(4))
    set(grob.hf,'position',[4 screensize(4)-pos.hf(4)-52 pos.hf(3:4)])
end

% Finish if we don't want to (re)define the positions
if ~resetflag && ~newgrob, return, end

% Make a vector of all objects - remember names and numbers
F = fieldnames(grob);
obj = [];
field = {};
num = [];
for k=1:length(F)
    f = F{k};
    nk = length(grob.(f));
    obj = [obj grob.(f)(:)'];
    field = [field repmat({f},1,nk)];
    num = [num 1:nk];
end
nobj = length(obj);

% Store properties which will be changed
state.units = fn_get(obj,'units');
state.buttondownfcn = fn_get(obj,'buttondownfcn');
hc=cell2mat(fn_get(obj,'children')); hc = hc(:)';
hctype = fn_get(hc,'type');
hc(fn_ismemberstr(hctype,{'uimenu','uicontextmenu'})) = [];
state.children = setdiff(hc,obj);
state.allobj = [obj state.children];
state.visible = fn_get(state.allobj,'visible');
type = fn_get(obj,'type');
state.figures = obj(fn_ismemberstr(type,'figure'));
state.figureresize     = fn_get(state.figures,'resize');
state.figureresizefcn  = fn_get(state.figures,'resizefcn');
state.figurewinbtdwnfcn = fn_get(state.figures,'windowbuttondownfcn');
state.axes = obj(fn_ismemberstr(type,'axes'));
state.axesaspect = fn_get(state.axes,'DataAspectRatio');
state.axesaspectmode = fn_get(state.axes,'DataAspectRatioMode');
state.uicontrols = obj(fn_ismemberstr(type,'uicontrol'));
state.uicontrolenable = fn_get(state.uicontrols,'enable');

% Scales
info.hf = grob.hf;
info.scales = [1 2 4 8 16 32 64 128 256];
scalek = 4; nscales = length(info.scales);

% Special buttons
pospanel = [2 102 96 24];
panel = uipanel('parent',info.hf,'units','pixels','position',pospanel, ...
    'buttondownfcn',@(hp,evnt)movepanel(hp,pospanel));
info.scale = pointer(info.scales(scalek));
info.scaleht = uicontrol('parent',panel,'style','text', ...
    'units','pixels','position',[21 1 30 16], ...
    'string',num2str(getvalue(info.scale)), ...
    'buttondownfcn',@(hu,evnt)movepanel(panel,pospanel),'enable','inactive');
uicontrol('parent',panel,'style','slider', ...
    'units','pixels','position',[1 1 20 20], ...
    'value',scalek,'sliderstep',[1 1]/(nscales-1),'min',1,'max',length(info.scales), ...
    'callback',@(hu,evnt)scaleupdate(info,hu)); % update value control
donehu = uicontrol('parent',panel,'style','pushbutton', ...
    'units','pixels','position',[53 1 40 20], ...
    'string','done','callback','delete(gcbo)');

% Change some properties
set(obj,'units','pixels','visible','on', ...
    'buttondownfcn',@(hobj,evnt)frameresize(info,hobj))
set(state.children,'visible','off')
set(state.figures,'buttondownfcn','','windowbuttondownfcn','', ...
    'resize','on','resizefcn',@(hf,evnt)figureresize(info,hf))
% set(state.axes,'dataaspectratiomode','auto')
for ha=state.axes, axis(ha,'normal'), end % strange that line above is not enough...
set(state.uicontrols,'enable','inactive')

% Wait
waitfor(donehu)

% Remove special buttons
delete(panel)

% Restore old property values
for k=1:nobj
    set(obj(k),'units',state.units{k},'buttondownfcn',state.buttondownfcn{k})
end
for k=1:length(state.allobj)
    set(state.allobj(k),'visible',state.visible{k})
end
for k=1:length(state.figures)
    set(state.figures(k),'windowbuttondownfcn',state.figurewinbtdwnfcn{k})
    set(state.figures(k),'resizefcn',state.figureresizefcn{k}), drawnow
    pos = get(state.figures(k),'pos'); % bug! changing 'resize' property changes the size!
    set(state.figures(k),'resize',state.figureresize{k}), drawnow
    set(state.figures(k),'pos',pos); 
end
for k=1:length(state.axes)
    set(state.axes(k),'dataaspectratio',state.axesaspect{k}, ...
        'dataaspectratiomode',state.axesaspectmode{k})
end
for k=1:length(state.uicontrols)
    set(state.uicontrols(k),'enable',state.uicontrolenable{k})
end

% Get position information
pos = [];
iobj = 0;
for k=1:length(F)
    f = F{k};
    nk = numel(grob.(f));
    positions = get(obj(iobj+(1:nk)),'position');
    if ~iscell(positions), positions = {positions}; end
    pos.(f) = cat(1,positions{:});
    iobj = iobj+nk;
end

% Command to generate pos
if nargout==0 || nargout==2
    cmd = [];
    for k=1:length(F)
        f = F{k};
        nk = length(grob.(f));
        cmdk = ['pos.' f ' = ['];
        for i=1:nk
            cmdk = [cmdk num2str(pos.(f)(i,:))]; %#ok<*AGROW>
            if i<nk, cmdk = [cmdk '; ']; end %#ok<*AGROW>
        end
        cmdk = [cmdk '];'];
        cmd = strvcat(cmd,cmdk); %#ok<VCAT>
    end
    if nargout==0
        disp(cmd)
        clear pos
    end
end


%-----------%
% CALLBACKS %
%-----------%

function scaleupdate(info,hu)

setvalue(info.scale,info.scales(get(hu,'value')));
set(info.scaleht,'string',num2str(getvalue(info.scale)))

%---
function figureresize(info,hf)

pos = get(info.hf,'position');
pos(3:4) = getvalue(info.scale)*(round(pos(3:4)/getvalue(info.scale)));
set(hf,'position',pos)

%---
function frameresize(info,hobj)

pos = get(hobj,'position');
hf = get(hobj,'parent');
posf = get(hf,'position');
p = get(0,'pointerlocation')-posf(1:2);
TOL = 2;
switch get(hf,'SelectionType')
    case 'normal'
        p = p-pos(1:2)+1.5;
        if p(2)<TOL
            if p(1)<TOL
                cat = 'botl';
                idx = {1 2};
            elseif p(1)>pos(3)-TOL
                cat = 'botr';
                idx = {3 2};
            else
                cat = 'bottom';
                idx = {[] 2};
            end
        elseif p(2)>pos(4)-TOL
            if p(1)<TOL
                cat = 'topl';
                idx = {1 4};
            elseif p(1)>pos(3)-TOL
                cat = 'topr';
                idx = {3 4};
            else
                cat = 'top';
                idx = {[] 4};
            end
        else
            if p(1)<TOL
                cat = 'left';
                idx = {1 []};
            elseif p(1)>pos(3)-TOL
                cat = 'right';
                idx = {3 []};
            else
                disp('please click a border of object')
                return
            end
        end
        set(info.hf,'pointer',cat)
        fn_buttonmotion({@frameresizeborder,hf,hobj,idx,getvalue(info.scale)},hf)
        set(hf,'pointer','arrow')
    case 'extend'
        p0 = getvalue(info.scale)*round(p/getvalue(info.scale));
        set(info.hf,'pointer','hand')
        fn_buttonmotion({@frameresizemove,hf,hobj,p0,pos,getvalue(info.scale)},hf)
        set(hf,'pointer','arrow')
    case 'alt'
        pos = getvalue(info.scale)*(round(pos/getvalue(info.scale)));
        set(hobj,'position',pos);
end

%---
function frameresizeborder(hf,hobj,idx,scale)

pos = get(hobj,'position');
pos(3:4) = pos(3:4)+pos(1:2);
posf = get(hf,'position');
p = get(0,'pointerlocation')-posf(1:2);
p = scale*round(p/scale);
pos(idx{1}) = p(1);
pos(idx{2}) = p(2);
pos(3:4) = max(1,pos(3:4)-pos(1:2));
set(hobj,'position',pos)

%---
function frameresizemove(hf,hobj,p0,pos,scale)

posf = get(hf,'position');
p = get(0,'pointerlocation')-posf(1:2);
p = scale*round(p/scale);
pos(1:2) = pos(1:2)+(p-p0);
set(hobj,'position',pos)

%---
function movepanel(hp,pos0)

hf = get(hp,'parent');
switch get(hf,'SelectionType')
    case 'alt'
        set(hp,'position',pos0)
    otherwise
        set(hf,'pointer','hand')
        fn_moveobject(hp)
        set(hf,'pointer','arrow')
end





