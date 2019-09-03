classdef fn_slider < hgsetget
    % function fn_slider([hp,][properties])
    %
    % See also fn_sliderenhance, fn_sensor, fn_control
    
    % Thomas Deneux
    % Copyright 2007-2012
    properties
        inc
        minmax = [0 1];
    end
    properties (Dependent)
        mode        % 'point', 'area', or 'area+point'
        value       % scalar in 'point' mode, 2-element vector in 'area' or 'area+point' mode
        point       % scalar in 'area+point' mode, [] otherwise
        width
    end
    properties (Dependent, GetAccess='private')
        sliderstep
        min
        max
    end
    properties
        callback
    end
    properties (Dependent)
        units
        position
        foregroundcolor
        backgroundcolor
        slidercolor
        visible
        deletefcn
    end
    properties (SetAccess='private')
        sliderscrolling = false;
        posframepix;
    end
    properties (Access='private')
        hf              % figure
        hpanel          % uipanel
        hframe          % frame around panel
        hslider         % uicontrol, frame style
        hline           % line, this is a second slider in fact
    end
    properties (Access='private')
        area = 1;
        x = [.25 .75];  % relative values; [left right] in 'area' mode, [value width] in 'point' mode
    end
    properties (Dependent, Access='private')
        sides
        left           
        right           
        center      
        pointpos
    end
    
    % Constructor/Destructor
    methods
        function U = fn_slider(varargin)
            % Objects
            if nargin>0
                a = varargin{1};
                if isscalar(a) && ishandle(a) && strcmp(get(a,'type'),'uipanel')
                    U.hpanel = a;
                    varargin(1)=[];
                end
            end
            if isempty(U.hpanel)
                % 'parent' property specified in arguments?
                hpar = [];
                for k=1:nargin
                    if isequal(varargin{k},'parent')
                        hpar = varargin{k+1};
                        varargin([k k+1])=[];
                        break
                    end
                end
                if isempty(hpar), hpar = gcf; end
                U.hpanel = uipanel('parent',hpar, ...
                    'units','pixel','position',[20 20 200 20]);
            end
            getnewframepos(U)
            set(U.hpanel,'bordertype','line','borderwidth',0, ...
                'resizefcn',@(u,evnt)getnewframepos(U))
            U.hframe = uicontrol('parent',U.hpanel, ...
                'style','frame','enable','off', ...
                'units','normalized','position',[0 0 1 1]);
            U.hslider = uicontrol('parent',U.hpanel, ...
                'style','frame','enable','off', ...
                'units','normalized');
            U.hline = uicontrol('parent',U.hpanel, ...
                'style','frame','enable','off','visible','off', ...
                'units','pixels');
            U.foregroundcolor = [0 0 0];
            U.slidercolor = [.3 .4 .5];
            U.hf = get(U.hpanel,'parent');
            while ~strcmp(get(U.hf,'type'),'figure'), U.hf = get(U.hf,'parent'); end
           
            % Position
            sliderposition(U)
            
            % Callbacks
            set(U.hframe,'buttondownfcn',@(u,evnt)event(U,'frame'))
            set(U.hslider,'buttondownfcn',@(u,evnt)event(U,'slider'))
            set(U.hline,'buttondownfcn',@(u,evnt)event(U,'line'))
            
            % User settings
            if ~isempty(varargin), set(U,varargin{:}), end
        end
        function delete(U)
            if ishandle(U.hpanel), delete(U.hpanel), end
        end
    end
    
    % GET/SET - basic object properties
    methods
        function str = get.mode(U)
            switch U.area
                case 0
                    str = 'point';
                case 1
                    str = 'area';
                case 2
                    str = 'area+point';
            end
        end
        function x = get.visible(U)
            x = get(U.hpanel,'visible');
        end
        function set.visible(U,x)
            set(U.hpanel,'visible',x)
        end
        function x = get.deletefcn(U)
            x = get(U.hpanel,'deletefcn');
        end
        function set.deletefcn(U,x)
            set(U.hpanel,'deletefcn',x)
        end
        function c = get.foregroundcolor(U)
            c = get(U.hframe,'foregroundcolor');
        end
        function set.foregroundcolor(U,c)
            set([U.hframe U.hslider U.hline],'foregroundcolor',c);
        end
        function c = get.backgroundcolor(U)
            c = get(U.hpanel,'backgroundcolor');
        end
        function set.backgroundcolor(U,c)
            set(U.hpanel,'backgroundcolor',c);
        end
        function c = get.slidercolor(U)
            c = get(U.hslider,'backgroundcolor');
        end
        function set.slidercolor(U,c)
            set(U.hslider,'backgroundcolor',c);
        end
        function pos = get.units(U)
            pos = get(U.hpanel,'units');
        end
        function set.units(U,pos)
            set(U.hpanel,'units',pos);
        end
        function pos = get.position(U)
            pos = get(U.hpanel,'position');
        end 
        function set.position(U,pos)
            set(U.hpanel,'position',pos);
            sliderposition(U)
        end
    end
    
    % GET/SET - positions
    methods
        function getnewframepos(U)
            tmp = get(U.hpanel,'units');
            set(U.hpanel,'units','pixel');
            U.posframepix = get(U.hpanel,'position');
            set(U.hpanel,'units',tmp)
        end
        function sliderposition(U)
            sid = U.sides;
            xpos = [sid(1) diff(sid)];
            if U.area
                if xpos(2)==0, xpos(2)=1e-4; end
            else
                % slider must have a minimal width of 10 pixels
                W = 10;
                if isvertical(U)
                    fact = U.posframepix(4);
                else
                    fact = U.posframepix(3);
                end
                w = W/fact;
                if xpos(2)<w;
                    xpos(1) = max(0,xpos(1)-w/2);
                    xpos(2) = w - max(0,w/2-xpos(1)) - max(0,w/2-(1-xpos(1))); 
                end
            end
            % actual position
            if isvertical(U)
                pos = [0 xpos(1) 1 xpos(2)];
            else
                pos = [xpos(1) 0 xpos(2) 1];
            end
            % update display if possible
            if diff(U.minmax)
                set(U.hslider,'position',pos)
            end
            % update line as well if mode area+point
            if U.area==2, lineposition(U), end
        end
        function lineposition(U)
            if U.area~=2, set(U.hline,'visible','off'), return, end
            xpos = U.pointpos;
            % line has a width of 2 pixels, so it is more convenient to
            % work in pixel units
            if isvertical(U)
                xpospix = xpos*U.posframepix(4);
                pos = [1 xpospix-1 U.posframepix(3) 2];
            else
                xpospix = xpos*U.posframepix(3);
                pos = [xpospix-1 1 2 U.posframepix(4)];
            end
            % update display if possible
            if diff(U.minmax)
                set(U.hline,'visible','on','position',pos)
            end
        end
        function set.x(U,x)
            % coerce
            if U.inc, step = 1/U.inc; end %#ok<*MCSUP>
            switch U.area
                case 0
                    % x(1) is relative value, x(2) is relative width
                    x(1) = max(0,min(1,x(1)));
                    if U.inc
                        x(1) = round(x(1)*step)/step;
                    end
                    if x(2)<=0 || x(2)>=1
                        error('relative width must be >0 and <1')
                    end
                case 1
                    % x(1) is left, x(2) is right
                    if diff(x)<0, x = x([2 1]); end
                    x = max(0,min(1,x));
                    if U.inc
                        x = round(x*step)/step;
                    end
                case 2
                    % x(1) is left, x(2) is right, x(3) is relative value
                    if diff(x(1:2))<0, x = x([2 1]); end
                    x = min(1,max(0,x));
                    % TODO: change set.mode so that the case that
                    % length(x)==2 in 'area+point' mode will not happen
                    if length(x)<3
                        x(3) = mean(x);
                    else
                        x(3) = max(x(1),min(x(2),x(3)));
                    end
                    if U.inc
                        x = round(x*step)/step;
                    end
            end
            U.x = x;
            % update display
            sliderposition(U)
        end
        function left = get.left(U)
            if U.area
                left = U.x(1);
            else
                left = U.x(1) * (1-U.x(2));
            end
        end
        function set.left(U,left)
            if U.area
                U.x(1) = left;
            else
                U.x(1) = left / (1-U.x(2));
            end
        end
        function right = get.right(U)
            if U.area
                right = U.x(2);
            else
                right = U.left + U.x(2);
            end
        end
        function set.right(U,right)
            if U.area
                U.x(2) = right;
            else
                error('cannot set right in ''point'' mode')
            end
        end
        function sides = get.sides(U)
            if U.area
                sides = U.x(1:2);
            else
                sides = U.left + [0 U.x(2)];
            end
        end
        function set.sides(U,sides)
            if U.area
                U.x(1:2) = sides;
            else 
                w = diff(sides);
                U.x = [sides(1)/(1-w) w];
            end
        end
        function c = get.center(U)
            if U.area
                c = mean(U.x);
            else
                c = U.left + U.x(2)/2;
            end
        end
        function set.center(U,c)
            if U.area
                w = diff(U.x(1:2));
                xnew = c + [-w/2 w/2];
                % do not allow width to change! (refuse to move center more
                % than a certain quantity)
                if xnew(1)<0
                    xnew = xnew-xnew(1);
                elseif xnew(2)>1
                    xnew = xnew-(xnew(2)-1);
                end
                U.x(1:2) = xnew;
            else
                U.left = c - U.x(2)/2;
            end
        end
        function p = get.pointpos(U)
            if U.area==2
                p = U.x(3);
            else
                p = [];
            end                
        end
        function set.pointpos(U,p)
            if U.area==2
                U.x(3) = p;
            else
                error('''pointpos'' property can be set only in ''area+point'' mode')
            end                
        end
    end
    
    % GET/SET - active properties
    methods
        function set.mode(U,str)
            if strcmp(str,U.mode), return, end
            % memorize position
            sid = U.sides;
            % update mode
            switch str
                case 'point'
                    U.area = 0;
                case 'area'
                    U.area = 1;
                case 'area+point'
                    U.area = 2;
                otherwise
                    error('wrong mode ''%s''',str)
            end
            % update value storage (U.x) by re-setting the position
            % this works in the case of 'area+point' mode, but it is
            % tricky, see set.x function
            U.sides = sid;
            % hide the line if needed
            if U.area~=2, set(U.hline,'visible','off'), end
        end     
        function set.inc(U,inc)
            if inc && mod(1,inc)
                error('rounding increment must divide 1')
            end
            U.inc = inc;
            U.x = U.x; % automatic update! (coerce)
        end
        function width = get.width(U)
            if U.area
                width = diff(U.x(1:2));
            else
                % moving the cursor by w is changing left by w, therefore
                % changes the value by w/(1-w)
                w = U.x(2);
                width = w/(1-w);
            end
        end
        function set.width(U,width)
            if U.inc
                % amount by which increase/decrease should be a multiple of
                % the rounding
                width = round(width/U.inc)*U.inc;
            end
            if U.area
                w = width;
                % try to preserve center
                c = U.center;
                xnew = c + [-w/2 w/2];
                if xnew(1)<0
                    xnew = xnew-xnew(1);
                elseif xnew(2)>1
                    xnew = xnew-(xnew(2)-1);
                end
                U.x(1:2) = xnew;
            else
                % if there are 1/width steps from min to max, there are
                % (1/width+1) positions in total, and the cursor relative
                % width should be 1/(1/width+1)
                w = width/(1+width);
                U.x(2) = w;
            end
        end
        function x = get.sliderstep(U)
            x = [U.inc U.width];
        end
        function set.sliderstep(U,x)
            U.inc = x(1);
            U.width = x(2);
        end
        function val = get.value(U)
            if U.area
                v = U.x(1:2);
            else
                v = U.x(1);
            end
            val = U.min + v*diff(U.minmax);
            % eliminate floating point error
            val = fn_round(val);
        end
        function set.value(U,val)
            if all(val==U.value), return, end
            v = (val-U.min) / diff(U.minmax);
            if U.area
                U.x = v;
            else
                U.x(1) = v;
            end
        end
        function val = get.point(U)
            if U.area==2
                v = U.pointpos;
            else
                v = [];
            end
            val = U.min + v*diff(U.minmax);
        end
        function set.point(U,val)
            if U.area~=2
                error('''pointpos'' property can be set only in ''area+point'' mode')
            end
            if all(val==U.point), return, end
            v = (val-U.min) / diff(U.minmax);
            U.pointpos = v;
        end
        function set.minmax(U,mM)
            if all(mM==U.minmax), return, end
            % memorize value
            val = U.value;
            p = U.point;
            % change min-max
            U.minmax = mM;
            % re-set value
            U.value = val;
            if U.area==2, U.point = p; end
        end
        function m = get.min(U)
            m = U.minmax(1);
        end
        function set.min(U,m)
            U.minmax(1) = m;
        end
        function M = get.max(U)
            M = U.minmax(2);
        end
        function set.max(U,M)
            U.minmax(2) = M;
        end
    end
    
    % Callbacks
    methods
        function event(U,flag)
            if strcmp(flag,'frame')
                % step the slider
                xf = mouseposframe(U);
                dir = 2*(xf>U.center)-1;
                U.center = U.center + dir*U.width;
            else 
                % slide
                p0 = mouseposframe(U);
                if strcmp(flag,'slider') && U.area
                    PIX = 2;
                    xs = mouseposslider(U);
                    if xs(1)<=PIX
                        flag = 'left';
                    elseif xs(2)>=-PIX
                        flag = 'right';
                    elseif U.area==2
                        % step the line
                        dir = 2*(p0>U.pointpos)-1;
                        step = fn_switch(U.inc~=0,U.inc,U.width/10);
                        U.pointpos = U.pointpos + dir*step;
                    end
                end
                switch flag
                    case 'left'
                        % move the left side of the slider
                        pos0 = U.left;
                        set(U.hf,'pointer',fn_switch(isvertical(U),'bottom','left'))
                    case 'right'
                        % move the right side of the slider
                        pos0 = U.right;
                        set(U.hf,'pointer',fn_switch(isvertical(U),'top','right'))
                    case 'slider'
                        % move the slider
                        pos0 = U.center;
                    case 'line'
                        % move the line
                        pos0 = U.pointpos;
                end
                U.sliderscrolling = true;
                fn_buttonmotion({@slide,U,flag,p0,pos0},U.hf)
                U.sliderscrolling = false;
                set(U.hf,'pointer','arrow')
            end
            exec(U)
        end
        function slide(U,flag,p0,pos0)
            p = mouseposframe(U);           
            switch flag
                case 'left'
                    U.left = p;
                case 'right'
                    U.right = p;
                case 'slider'
                    U.center = pos0+(p-p0);
                case 'line'
                    U.pointpos = pos0+(p-p0);
            end
            exec(U)
        end
        function exec(U)
            if isempty(U.callback), return, end
            switch class(U.callback)
                case 'char'
                    evalin('base',U.callback)
                case 'function_handle'
                    feval(U.callback,U,[]);
                case 'cell'
                    feval(U.callback{1},U,[],U.callback{2:end})
            end
        end
    end
    
    % Routines
    methods
        function b = isvertical(U)
            pos = U.posframepix;
            b = pos(4)>pos(3);
        end
        function x = mouseposframe(U)
            % normalized position inside frame
            pos = mousepos(U.hpanel);
            if isvertical(U)
                x = (pos(2)-1) / U.posframepix(4);
            else
                x = (pos(1)-1) / U.posframepix(3);
            end
        end
        function x = mouseposslider(U)
            % pixel positions relative to left and right sides of slider
            pos = mousepos(U.hpanel);
            set(U.hslider,'units','pixels')
            posframe = get(U.hslider,'position');
            set(U.hslider,'units','normalized')
            if isvertical(U)
                x = [pos(2)-posframe(2) pos(2)-(posframe(2)+posframe(4))];
                x = x+1; % BERK
            else
                x = [pos(1)-posframe(1) pos(1)-(posframe(1)+posframe(3))];
            end
        end
    end
    
end

function pos = mousepos(hobj)
% position in pixel units of pointer in current container (figure or
% uipanel)
% this is SHITTY!!!

switch get(hobj,'type')
    case 'figure'
        tmp = get(hobj,'units');
        set(hobj,'units','pixel')
        pos = get(hobj,'currentpoint');
        set(hobj,'units',tmp)
    case 'uipanel'
        tmp = get(hobj,'units');
        set(hobj,'units','pixel')
        panelpos = get(hobj,'position');
        set(hobj,'units',tmp)
        pos = mousepos(get(hobj,'parent'));
        pos = pos-(panelpos(1:2)-1);
    otherwise
        error programming
end

end