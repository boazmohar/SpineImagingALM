function mask = fn_maskselect(a,mouseflag,dorepeat)
% function mask = fn_maskselect(image[,mouseflag][,dorepeat])
%---
% 
% Input:
% - image       2D array
% - mouseflag   'rect', 'poly' [default], 'free', 'ellipse'
% - dorepeat    select multiple regions? [default = false]
% 
% Output:
% - mask        logical array the same size of image indicating interior of
%               the mask

% Thomas Deneux
% Copyright 2011-2012

% Input
[nx ny] = size(a);
if nargin<2, mouseflag = 'poly'; end
if nargin<3, dorepeat = false; end
switch mouseflag
    case 'rect'
        mouseflag = 'rectangle';
    case 'ellipse'
        error('not implemented yet')
    case {'poly' 'free'}
        % ok
    otherwise
        error('unknown flag ''%s''',mouseflag)
end
mouseflag= [mouseflag '+'];
    
% prepare display
hf = figure; set(hf,'tag','fn_maskselect')
imagesc(permute(a,[2 1 3])); axis image
if dorepeat
    more = uicontrol('string','more','style','togglebutton', ...
        'pos',[5 5 40 18]);
    uicontrol('string','ok','callback',@(u,e)close(hf), ...
        'pos',[50 5 20 18]);
end

% go (loop if dorepeat is true)
firsttime = true;
mask = false;
while true
    poly = fn_mouse(mouseflag);
    mask = mask | poly2mask(poly(2,:),poly(1,:),ny,nx);
    if dorepeat
        waitfor(more,'value',1)
        if ~ishandle(more), break, end
        set(more,'value',0)
    else
        close(hf)
        break
    end
end

