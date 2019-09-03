function m = fn_menu(h,varargin)
% function m = fn_menu(fig,'h'|'v',title,width,height,'FramePropertyName',FramePropertyValue,...)
% function u = fn_menu(m,'add','UIControlPropertyName',UIControlPropertyValue,...)
% function fn_menu(m,'redraw')
%---
% creates a frame uicontrol and provides facilities to put new uicontrols
% inside

% Thomas Deneux
% Copyright 2006-2012

if nargin<1
    h = gcf;
end
if mod(h,1)==0 && ~ishandle(h), figure(h), end

switch get(h,'type')

    case 'figure'       % NEW MENU
        if nargin<2, align='v'; else align=varargin{1}; end
        if nargin<3, title=''; else title=varargin{2}; end
        if nargin<4, width=30; else width=varargin{3}; end
        if nargin<5, height=20; else height=varargin{4}; end
        switch align
            case {'v','h'}
                position = [5 5 width+10 height+10];
            otherwise
                error('unknown align flag ''%s''',align)
        end
        m = uicontrol('parent',h,'style','frame','position',position,varargin{5:end},...
            'deletefcn',{@kill});
        updateinfo(m,struct('align',align,'width',width,'height',height,...
            'title',title,'children',[]));

    case 'uicontrol'    % EDIT MENU
        command = varargin{1};
        m = h;
        h = get(m,'parent'); % parent figure
        info = getappdata(m);
        switch command
            case 'add'
                u = uicontrol('parent',h,'visible','off',varargin{2:end});
                info.children(end+1) = u;
                % update positions
                set(m,'Position',controlposition(0,info))
                for i = 1:length(info.children)
                    set(info.children(i),'Position',controlposition(i,info))
                end
                set(u,'visible','on')
            case 'redraw'
                set(m,'Position',controlposition(0,info))
                for i = 1:length(info.children)
                    set(info.children(i),'Position',controlposition(i,info))
                end
            otherwise
                error('unknown command flag ''%s''',command)
        end
        updateinfo(m,info)
        clear m
        if nargout==1, m=u; end
end


% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%

function position = controlposition(i,info)

n = length(info.children);
switch info.align
    case 'v'    % vertical alignment
        if i==0     % frame
            position = [5 5 (info.width+10) (n*info.height+max(0,n-1)*2+10)];
        else        % child uicontrol
            position = [10 (10+(n-i)*(info.height+2)) info.width info.height];
        end
     case 'h'    % horizontal alignment
        if i==0     % frame
            position = [5 5 (n*info.width+max(0,n-1)*2+10) (info.height+10)];
        else        % child uicontrol
            position = [(10+(i-1)*(info.width+2)) 10 info.width info.height];
        end
end
       
%---
function updateinfo(m,info)

for f=fieldnames(info)'
    setappdata(m,f{1},info.(f{1}))
end

%---
function kill(m,dum)

info = getappdata(m);
delete(info.children(ishandle(info.children)))




