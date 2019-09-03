function fn_playmovie(M)
% function fn_playmovie(M)
%---
% Simple display of movie
% It is preferable to use fn_movie 
%
% See also fn_movie

% Thomas Deneux
% Copyright 2005-2012

n = size(M,3);
colormap gray

clf
im = imagesc(M(:,:,1)'); axis image
ha = axes('position',[.1 .05 .8 .01],'box','on','xlim',[1 n],'ylim',[-1 1],'ytick',[]);
pt = line('parent',ha,'xdata',1,'ydata',0,'marker','.','color','b','markersize',20, ...
    'buttondownfcn',{@pointhit});
set(ha,'buttondownfcn',{@axeshit,pt});
setappdata(pt,'mode','normal')
setappdata(pt,'play',true)

% play
while true
    if strcmp(getappdata(pt,'mode'),'exit')
        set(ha,'buttondownfcn','')
        set(pt,'buttondownfcn','')
        return
    end
    i = mod(get(pt,'xdata')-1,n)+1;
    if getappdata(pt,'play'), i = mod(i,n)+1; end
    set(pt,'xdata',i)
    set(im,'CData',M(:,:,i)')
    pause(.001)
end

%---
function pointhit(pt,dum)

if strcmp(get(gcf,'selectiontype'),'alt')
    setappdata(pt,'mode','exit')
    return
end

setappdata(pt,'mode','hit')
set(gcf,'windowbuttonmotionfcn',{@pointmotion,pt})
set(gcf,'windowbuttonupfcn',{@pointrelease,pt})

%---
function pointmotion(dum1,dum2,pt)

setappdata(pt,'mode','motion')

p = get(gca,'currentpoint'); i = round(p(1));
set(pt,'xdata',i)

%---
function pointrelease(dum1,dum2,pt)

set(gcf,'windowbuttonmotionfcn','')
set(gcf,'windowbuttonupfcn','')

play = getappdata(pt,'play');
if strcmp(getappdata(pt,'mode'),'hit')  
    % point has not been moved -> pause toggle
    play = ~play;
end
setappdata(pt,'play',play)
if play
    set(pt,'color','b')
else
    set(pt,'color','r')
end

setappdata(pt,'mode','normal')

%---
function axeshit(dum1,dum2,pt)

if strcmp(get(gcf,'selectiontype'),'alt')
    setappdata(pt,'mode','exit')
    return
end

if getappdata(pt,'play')
    setappdata(pt,'play',false)
    set(pt,'color','r')
else
    i = get(pt,'xdata');
    p = get(gca,'currentpoint'); x = p(1);
    if x<i, i=max(1,i-1); else xlim = get(gca,'xlim'); i=min(xlim(2),i+1); end
    set(pt,'xdata',i)
end


