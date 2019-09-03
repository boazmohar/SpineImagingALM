function fn_showcolormap(cmap,varargin)
% function fn_showcolormap(cmap[,clip][,ha])
%---
% opens a new figure and display in it the colormap at the given clipping
% value

% Thomas Deneux
% Copyright 2011-2012

if nargin==0, help fn_showcolormap, return, end

% input
if ischar(cmap), cmap = feval(cmap,256); end
defclip = false; ha = [];
for k=1:length(varargin)
    a = varargin{k};
    if ishandle(a)
        ha = a;
    else
        defclip = true;
        clip = a;
    end
end
if ~defclip, clip = [0 1]; end
if isempty(ha)
    hf = figure(873);
    clf(hf), fn_figmenu(hf)
    fn_setfigsize(hf,200,500);
    set(hf,'color','w','tag','color bar')
    ha = axes('pos',[.5 .2 .1 .7]);
end

% display
colormap(ha,cmap)
imagesc([0 1],clip,linspace(clip(1),clip(end),size(cmap,1))','parent',ha)
axis(ha,'normal')
set(ha,'ydir','normal','xtick',[])
if ~defclip, set(ha,'ytick',[]), end

