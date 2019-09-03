function fn_controlpositions(hu,ha,posa,pospix)
% function fn_controlpositions(hu,ha,posa,pospix)
%---
% set the position of controls relatively to an axes, and set the figure
% resize function to update these positions automatically
%
% Input
% - hu      control handle
% - ha      axes handle
% - posa    position relative to axes ([0 0] = bottom-left corner, [1 1] =
%           up-right corner)
% - pospix  position in pixels to add to 'posa' and size of control

% Thomas Deneux
% Copyright 2007-2012

% attach information to controls
info = struct('ha',ha,'posa',posa,'pospix',pospix);
setappdata(hu,'fn_controlpositions',info);

% set figure resize function
hf = get(ha,'parent');
fun = get(hf,'resizefcn');
if ~isempty(fun) && ~isequal(fun,@updatepositions)
    error('a resize function already exists in figure')
end
set(hf,'resizefcn',@updatepositions,'doublebuffer','on')

% update for this control
updatepositions(hf,[],hu)

%---
function updatepositions(hf,evnt,hlist) 

if nargin<3
    hlist = findobj(hf,'type','uicontrol');
end

for hu = hlist'    
    info = getappdata(hu,'fn_controlpositions');
    if isempty(info), continue, end
    posa = fn_coordinates(info.ha,'c2g',info.posa(1:2),'position')';
    if length(info.posa)==4
        posa(3:4) = fn_coordinates(info.ha,'c2g',info.posa(3:4),'vector')';
    else
        posa(3:4) = 0;
    end
    pos = posa+info.pospix;
    pos([3 4]) = max(pos([3 4]),2);
    set(hu,'position',pos)
end

