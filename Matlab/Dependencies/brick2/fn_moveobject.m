function dp = fn_moveobject(hobj,varargin)
% function dp = fn_moveobject(hobj[,'fast'|('fastcolor',col)][,'latch'][,'point',i])
%---
% moves objects while mouse button is pressed
% 
% Options
% - 'fast'      use 'xor' erase mode for faster display updates
% - 'fastcolor' color to use in 'xor' mode
% - 'latch'     when button is released, brings objects back to initial
%               position
% - 'point',i   move only the ith point(s) of line objects
% - 'twice'     wait for button press+release, or release+pressagain
% 
% See also fn_buttonmotion

% Thomas Deneux
% Copyright 2007-2012

% Input
if ~isvector(hobj) || ~all(ishandle(hobj))
    error('hobj must be a vector of graphic handles')
end
par = fn_get(hobj,'parent');
if ~isscalar(unique([par{:}])), error('objects must have the same parent'), end
% (options)
fastflag = false; latchflag = false; twiceflag = false; pointidx = 0;
i = 1; 
while i<=length(varargin)
    switch varargin{i}
        case 'fast'
            fastflag = true;
            fastcolor = [1 1 1]*.5;
        case 'fastcolor'
            fastflag = true;
            fastcolor = varargin{i+1}; i=i+1;            
        case 'latch'
            latchflag = true;
        case 'point'
            pointidx = varargin{i+1}; i=i+1;
        case 'twice'
            twiceflag = true;
        otherwise
            error('unknown flag''%''',varargin{i})
    end
    i = i+1;
end

% Referential, units, property changes
switch get(get(hobj(1),'parent'),'type')
    case {'figure','uipanel'}
        % pointer location will be in screen referential (pixel units), so
        % change 'units' property of objects to 'pixel'
        ref = 0;
        oldunits = fn_get(hobj,'units');
        set(hobj,'units','pixels')
    case {'axes'}
        % pointer location will be an axes referential
        ref = get(hobj(1),'parent');
        if fastflag
            oldstate = fn_get(hobj,{'color','erasemode'});
            fn_set(hobj,{'color','erasemode'},{fastcolor,'xor'})
        end
    otherwise
        error('cannot handle parent type ''%s''',get(get(hobj(1),'parent'),'type'))
end

% Find parent figure
hf = get(hobj(1),'parent');
while ~strcmp(get(hf,'type'),'figure'), hf = get(hf,'parent'); end

% Moving object while button is pressed
p0 = getpoint(ref);
pos0 = getpos(hobj,ref);
fn_buttonmotion({@movesub,hobj,ref,p0,pos0,pointidx},hf)
if twiceflag
    fn_buttonmotion({@movesub,hobj,ref,p0,pos0,pointidx},hf)
end

% Restore properties
if latchflag, setpos(hobj,ref,p0,pos0,p0,pointidx), end
if ~ref
    fn_set(hobj,'units',oldunits)
else
    if fastflag
        fn_set(hobj,oldstate)
    end
end

% output
if nargout>0
    p = getpoint(ref);
    dp = p-p0;
end

%---
function movesub(hobj,ref,p0,pos0,pointidx)

p = getpoint(ref);
setpos(hobj,ref,p0,pos0,p,pointidx)

%---
function p = getpoint(ref)

if ~ref
    p = get(ref,'pointerlocation');
else
    p = get(ref,'currentpoint');
    p = p(1,1:2);
end

%---
function pos = getpos(hobj,ref)

if ~ref
    pos = fn_get(hobj,'position');
else
    nobj = length(hobj);
    pos = cell(1,nobj);
    for k=1:nobj
        if isprop(hobj(k),'position')
            pos{k} = get(hobj(k),'position');
        elseif isprop(hobj(k),'xdata')
            pos{k} = get(hobj(k),{'xdata','ydata'});
        else
            error('don''t know how to handle ''%s'' object', ...
                get(hobj(k),'type'))
        end
    end
end

%---
function setpos(hobj,ref,p0,pos0,p,pointidx)

dp = p(1,1:2)-p0;
if ~ref
    for k=1:length(hobj)
        posk = pos0{k};
        posk(1:2) = posk(1:2) + dp;
        set(hobj(k),'position',posk)
    end
else
    for k=1:length(hobj)
        posk = pos0{k};
        if isnumeric(posk)
            posk(1:2) = posk(1:2) + dp;
            set(hobj(k),'position',posk)
        else
            if pointidx
                posk{1}(pointidx)=posk{1}(pointidx)+dp(1);
                posk{2}(pointidx)=posk{2}(pointidx)+dp(2);
            else
                posk{1}=posk{1}+dp(1);
                posk{2}=posk{2}+dp(2);
            end
            set(hobj(k),{'xdata' 'ydata'},posk)
        end
    end
end


