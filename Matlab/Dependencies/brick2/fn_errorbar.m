function hl = fn_errorbar(varargin)
% function hl = fn_errorbar([x,]y,e,line options[,flag])
% function hl = fn_errorbar([x,]y,line options[,flag])
%---
% ym and ystd can be vectors or arrays
% if e is not supplied, y is redefined and e is defined as the mean and
% std/sqrt(n) of y along the 2d dimension (if y is a 2D array) or the 3d
% dimension (if y is a 3D array)
%
% flag can be 'bar' or 'patch'
%
% uses bar display instead of plot if 'bar' flag is specified

% Thomas Deneux
% Copyright 2006-2012

if nargin==0, help fn_errorbar, return, end

% Input
% ('bar' flag)
a = varargin{end};
if ischar(a) && fn_ismemberstr(a,{'bar' 'patch' 'xerror'})
    flag = a;
    varargin(end)=[];
    narg = nargin-1;
else
    flag = '';
    narg = nargin;
end
% (options)
for i=1:narg
    if i==narg || ischar(varargin{i+1}), break, end
end
narg = i;
opt = varargin(narg+1:end);
% (parent axes)
ha = [];
for i=1:length(opt)-1
    if ischar(opt{i}) && strcmp(opt{i},'parent')
        ha = opt{i+1};
    end
end
if isempty(ha), ha = gca; end
% (numeric arguments)
x = []; y = []; e = []; ex = [];
switch narg
    case 1
        Y = varargin{1};
        dim = ndims(Y);
        y   = mean(Y,dim);
        e = std(Y,0,dim)/sqrt(size(Y,dim));
    case 2
        if strcmp(flag,'xerror')
            [X Y] = deal(varargin{1:2});
            if ndims(X)~=2 || ndims(Y)~=2, error 'if error on x and y, >2d is not possible', end
            x  = mean(X,2);
            ex = std(X,0,2)/sqrt(size(X,2));
            y  = mean(Y,2);
            e  = std(Y,0,2)/sqrt(size(Y,2));
        elseif ~isvector(varargin{1}) || isvector(varargin{2})
            [y e] = deal(varargin{1:2});
        else
            [x Y] = deal(varargin{1:2});
            dim = ndims(Y);
            y   = mean(Y,dim);
            e = std(Y,0,dim)/sqrt(size(Y,dim));
        end
    case 3
        [x y e] = deal(varargin{1:3});
    case 4
        [x ex y e] = deal(varargin{1:4});
    otherwise
        error arguments
end
if isvector(y), y = y(:); e = e(:); end
if isempty(x), x = (1:size(y,1))'; else x = x(:); end

% Prepare for display
[nt n] = size(y);
cols = get(ha,'ColorOrder'); ncol = size(cols,1);
yb = y-e; yt = y+e;

% Display
switch flag
    case 'bar'
        % bar display
        ny = size(y,2);
        ddx = diff(x,2);
        if max(abs(ddx)) > 100*eps, error('x points not equidistant'); end
        dx = (x(2)-x(1)) / (ny+1.5);
        xdispatch = dx * (-(ny-1)/2 + (0:ny-1));
        xx = fn_add(x,xdispatch);
        isholdoff = ~strcmp(get(ha,'nextplot'),'add');
        if isholdoff
            errorbar(xx,y,e,'color','k','parent',ha);
            ax = axis(ha);
        end
        hl{1} = bar(x,y,opt{:});
        hold(ha,'on')
        htemp = errorbar(xx,y,e,'color','k','parent',ha);
        if isholdoff, hold(ha,'off'), end
        htemp = get(htemp,'children');
        if ny==1, htemp = {htemp}; end
        hl{2} = zeros(1,ny);
        for i=1:ny
            delete(htemp{i}(1))
            hl{2}(i) = htemp{i}(2); %#ok<AGROW>
        end
        if any(y(:)<0), m=ax(3); else m=0; end
        if isholdoff, axis(ha,[x(1)-dx*ny/2 x(end)+dx*ny/2 m ax(4)]), end
    case 'patch'
        % display error as thick bands
        hl = plot(x,[yb yt],'parent',ha); % first set the axis by calling plot
        delete(hl)
        hl = {zeros(1,n) zeros(1,n)};
        for k=1:n
            kc = 1+mod(k-1,ncol);
            hl{1}(k) = patch([x(1:nt)' x(nt:-1:1)'],[yb(:,k)' yt(nt:-1:1,k)'], ...
                (1+cols(kc,:))/2, ...
                'parent',ha, ...
                'edgecolor','none'); %,'facealpha',.5);
        end
        for k=1:n
            hl{2}(k) = line(x,y(:,k),'color',cols(k,:),opt{:},'parent',ha);
        end
    case 'xerror'
        % scatter plot with errors for both x and y
        nx = size(x,1);
        xdata = [x-ex x+ex nan(nx,1) x x nan(nx,1)]';
        ydata = [y y nan(nx,1) y-e y+e nan(nx,1)]';
        hl = zeros(2,nx);
        hl(2,:) = plot(xdata,ydata,'color','b','parent',ha);
        hl(1,:) = line(repmat(x',2,1),repmat(y',2,1),'color','b','parent',ha, ...
            'linestyle','none','marker','.','markersize',16);
    otherwise
        % display error with dotted lines
        hl = plot(x,[y yt yb],opt{:});
        % automatic colors and line style
        ncol = size(cols,1);
        for k=1:n
            set(hl(k+[0 n 2*n]),'color',cols(1+mod(k-1,ncol),:))
        end
        set(hl(1:n),'linestyle','-')
        set(hl(n+1:3*n),'linestyle','--')
        % user options
        nopt = length(opt);
        if nopt>=2
            % ignore first optional argument if the total number of arguments is
            % odd
            if mod(nopt,2)
                set(hl,opt{2:nopt})
            else
                set(hl,opt{:})
            end
        end
end

if nargout==0, clear hl, end    
