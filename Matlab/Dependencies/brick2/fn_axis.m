function axout = fn_axis(varargin)
% function ax = fn_axis([ha,]'tight|image|tightimage',sidefactor,'y0')
%---
% set a nice range to axes
%
% Input:
% - flag            'tight'     stretch the axis as much as possible
%                   'image'     set an equal ratio between x and y
%                   'tightimag' stretch the axis as much as possible, while
%                               maintaining an equal ratio between x and y
% - sidefactor      scalar or 2-elements vector >1 but close to one: the
%                   difference with 1 indicates the small gap to leave to 
%                   the side
% - 'y0' flag       set ymin to 0
% 
% Example:
%   plot(sin(0:.01:100))
%   fn_axis('tight',[1 1.2])

% Thomas Deneux
% Copyright 2010-2012

% Input
ha = []; flag = 'tight'; doy0 = false; sidefactor = [];
for k=1:nargin
    x = varargin{k};
    if ishandle(x)
        ha = x;
    elseif ischar(x)
        if strcmp(x,'y0'), doy0 = true; else flag = x; end
    else
        sidefactor = x;
    end
end
dotight = false; dosquare = false;
switch flag
    case 'tight'
        dotight  = true;
    case 'image'
        dosquare = true;
    case 'imagetight'
        dotight  = true;
        dosquare = true;
    otherwise
        error argument
end
if isempty(ha), ha = gca; end
if isempty(sidefactor), sidefactor = 1 + .2*dotight; end
if isscalar(sidefactor), sidefactor = [1 1]*sidefactor; end

% Make axis tight and get the range
if dotight, axis(ha,'tight'), end
ax = axis(ha);
if abs(ax(3))<abs(ax(4))*1e-6, doy0 = true; end
wh = [ax(2)-ax(1) ax(4)-ax(3)];
center = ax([1 3]) + wh/2;

% Modify it to satisfy the 'dosquare' flag if necessary, and leave the gap
% specified in sidefactor
factx = sidefactor(1); facty = sidefactor(2);
if dosquare
    oldunits = get(ha,'units');
    set(ha,'units','pixel')
    pos = get(ha,'pos'); 
    set(ha,'units',oldunits)
    r = (factx*pos(4))/(facty*pos(3));
    r1 = wh(2)/wh(1);
    if r1<r
        % increase y view
        facty = facty * (r/r1);
    else
        % increase x view
        factx = factx * (r1/r);
    end
end
wh = wh .* [factx facty];
if doy0
    % special case: set ymin to 0 
    ax = [center(1)+wh(1)/2*[-1 1] [0 center(2)+wh(2)/2]];
else
    ax = [center(1)+wh(1)/2*[-1 1] center(2)+wh(2)/2*[-1 1]];
end
axis(ha,ax)

if nargout>0, axout = ax; end
