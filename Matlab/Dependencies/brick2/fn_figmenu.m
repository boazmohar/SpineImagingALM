function fn_figmenu(hf,varargin)
% function fn_figmenu
% function fn_figmenu(hf)
%---
% adds a custom menu to figure
%
% Note that function fn_figmenu makes use of fn_cd, and in order to be
% able to use it to save figures in the predefined directory 'capture',
% you should use 'fn_cd edit' to define where is this directory 'capture'.

% Thomas Deneux
% Copyright 2007-2012

if nargin==0
    hfs = findobj('type','figure');
    for hf = hfs'
        m = findobj(hf,'tag','fn_figmenu');
        if isempty(m)
            m = uimenu(hf,'label','Utils','Tag','fn_figmenu'); 
        else
            delete(get(m,'children'))
        end
        setsubmenus(m,hf)
    end
else
    % check if figure is likely to be a dialog -> in such case, cancel the
    % menu creation
    if ~isempty(get(hf,'buttonDownFcn')), return, end
    m = uimenu(hf,'label','Utils','Tag','fn_figmenu');
    setsubmenus(m,hf)
end

%---
function setsubmenus(m,hf)

uimenu(m,'label','tmp','Callback','tmp')
uimenu(m,'label','edit tmp','Callback','edit tmp')
uimenu(m,'label','edit fn_figmenu','Callback','edit fn_figmenu')

uimenu(m,'label','fn_imvalue','Callback','fn_imvalue','separator','on')
uimenu(m,'label','fn_imvalue image','Callback','fn_imvalue image')
uimenu(m,'label','fn_imvalue clean','Callback','fn_imvalue clean')
uimenu(m,'label','fn_imvalue end','Callback','fn_imvalue end')
m1 = uimenu(m,'label','more');
uimenu(m1,'label','fn_imvalue xy image','Callback','fn_imvalue xy image')
uimenu(m1,'label','fn_imvalue register ...','Callback','fn_imvalue register')
uimenu(m1,'label','fn_imvalue unregister','Callback','fn_imvalue unregister')

uimenu(m,'label','colormap gray','Callback','colormap gray','separator','on')
uimenu(m,'label','colormap jet','Callback','colormap jet')
m1 = uimenu(m,'label','more');
uimenu(m1,'label','colormap mapclip','Callback','colormap mapclip')
uimenu(m1,'label','colormap signcheck','Callback','colormap signcheck')
uimenu(m1,'label','colormap green','Callback','colormap green')

uimenu(m,'label','reset figure callbacks', ...
    'Callback', ...
    'set(gcf,''WindowButtonMotionFcn'','''',''WindowButtonUpFcn'','''',''KeyPressFcn'','''')', ...
    'separator','on')

uimenu(m,'label','fn_imdistline','Callback','fn_imdistline','separator','on')

uimenu(m,'label','save figure in PNG','Callback',@(u,evnt)savepng(hf,false,true),'separator','on')
uimenu(m,'label','save figure in PNG...','Callback',@(u,evnt)savepng(hf,false,false))
m1 = uimenu(m,'label','more');
uimenu(m1,'label','save image in PNG','Callback',@(u,evnt)savepng(hf,true,true))
uimenu(m1,'label','save image in PNG...','Callback',@(u,evnt)savepng(hf,true,false))

%---
function savepng(hf,doimg,doautofile)

base = fn_switch(doimg,'image','figure');
if doautofile
    fname = [fn_cd('capture') '/' base num2str(hf) '_' datestr(now,'YYmmDDHHMMSS') '.png'];
else
    fname = fn_savefile;
    [b fil ext] = fileparts(fname);
    if isempty(ext)
        fname = [b '/' fil '.png'];
    end
end
if doimg
    ha = get(hf,'children');
    if ~isscalar(ha), ha = gca; end
    c = get(ha,'children');
    f = find(strcmp(get(c,'type'),'image'));
    if isempty(f), return, end
    a = get(c(f(1)),'cdata');
    if ndims(a)==2
        cmap = get(hf,'colormap');
        clip = get(ha,'clim');
        fn_saveimg(a',fname,clip,1,cmap)
    else
        fn_saveimg(permute(a,[2 1 3]),fname)
    end
else
    fn_savefig(hf,fname)
end


