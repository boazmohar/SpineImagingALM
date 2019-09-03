function shiftman = fn_alignimage(a,b,varargin)
% function [shift =] fn_alignimage(a,b[,shift][,hf])
%---
% display a superposition of gray-level images a and b, and allow to move
% image b
% 
% Example:
%   load trees, 
%   Y = interp2(1:350,(1:258)',X,-4.3:344.7,(13:270)');
%   shift=fn_alignimage(X,Y)

% Thomas Deneux
% Copyright 2011-2012

% Input
a = double(a);
b = double(b);
[nx ny] = size(a);
if any(size(b)~=[nx ny]), error 'dimension mismatch', end
shift = [0 0]; hf = [];
for k=1:length(varargin)
    x = varargin{k};
    if length(x)>1
        shift = x;
    else
        hf = x;
    end
end
if isempty(hf), hf = gcf; end

% Rescale
a = a-min(a(:)); a = a/max(a(:));
b = b-min(b(:)); b = b/max(b(:)); 

% Display
clf(hf)
set(hf,'tag','xx')
ha = axes('parent',hf);
hi = imagesc(cat(3,a,a,b),'parent',ha);
set(hi,'buttondownfcn',@(u,e)startmove,'hittest','on')
set(ha,'clim',[-.5 .5],'climmode','manual')
p0 = []; shift0 = [0 0];

hu1 = uicontrol('pos',[5 5 75 18],'parent',hf,'style','text', ...
    'callback',@(u,e)getshiftvalue);
hu = uicontrol('pos',[5 25 80 18],'parent',hf,'style','radiobutton','string','show diff', ...
    'callback',@(u,e)showimages);
showimages

set([hf hu1 hu],'keyPressFcn',@(u,e)keypress(e))

if nargout>=1
    figure(hf)
    fn_okbutton('wait','parent',hf)
    shiftman = shift;
end

    function startmove
        p0 = get(ha,'currentPoint'); p0 = p0(1,1:2);
        shift0 = shift;
        fn_buttonmotion(@moveb,hf)
    end
    function moveb
        p = get(ha,'currentPoint'); p = p(1,1:2);
        shift = shift0 + [p(2)-p0(2) p(1)-p0(1)];
        showimages
    end
    function keypress(e)
        d = [];
        switch e.Key
            case 'leftarrow'
                d = [0 -1];
            case 'rightarrow'
                d = [0 1];
            case 'downarrow'
                d = [1 0];
            case 'uparrow'
                d = [-1 0];
            case 'return'
                set(hu,'value',~get(hu,'value'))
                showimages
        end
        if isempty(d), return, end
        if ~isempty(e.Modifier)
            switch e.Modifier{1}
                case 'shift'
                    d = d*5;
                case 'control'
                    d = d/10;
            end
        end
        shift = shift+d;
        showimages
    end
    function showimages
        bmov = min(max(fn_translate(b,shift),0),1);
        if get(hu,'value')
            im = a-bmov;
        else
            im = cat(3,a,a,bmov);
        end
        set(hi,'cdata',im)
        set(hu1,'string',sprintf('x: %.2f, y: %.2f',shift(2),shift(1)));
    end

end

