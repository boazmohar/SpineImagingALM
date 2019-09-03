function y = fn_framedisplay(x,varargin)
% [y =] function fn_framedisplay(x[,clip][,'singleaxes|multaxes'][,hf|ha][,'display'][,ncol])
%---
% Input (except for x, order of inputs can be changed):
% - x       3D or 4D data; NOT using Matlab image convention (i.e. first
%           dim is x and second dim is y) 
% - clip    clipping range
% - axflag  'singleaxes' will create a large image by putting each frame
%           next to each other and display this image
%           'multaxes' will display each frame in a separate axes
% - hf|ha   handle of figure (implicit 'multaxes' option) or axes (implicit
%           'singleaxes' option)
% - 'display'   by default, nothing is displayed if an output is requested,
%           this option is used to display even if there is an output
% - ncol    specify the number of columns to use
%
% Output:
% - y       a large image obtained by concatenating all images, note that
%           compared to the large image created with the 'singleaxes'
%           option, this image is slighlty larger since there is a 1-pixel
%           wide separation between all frames

% Thomas Deneux
% Copyright 2011-2012


% Input
% (data)
if iscell(x), x = cat(3,x{:}); end
s = size(x); s(end+1:5) = 1;
nx = s(1); ny = s(2);
if isreal(x) && (s(4)==1 || s(3)~=3)
    % standard movie
    nc = 1;
    nt = [s(3) prod(s(4:end))];
elseif ~isreal(x)
    % complex numbers: real and imaginary part make 2 separate channels
    x = cat(3,real(x),imag(x));
    nc = 1;
    nt = [2 prod(s(3:end))];
elseif s(3)==3
    % color movie
    nc = 3;
    nt = [s(4) prod(s(5:end))];
else
    error('number of dimensions of x must be 3 or 4')
end
x = reshape(x,[nx ny nc nt]);
% (options)
clip = []; axesflag = ''; ha = []; 
dooutput = (nargout==1); dodisplay = ~dooutput;
ncol = [];
for k=1:length(varargin)
    a = varargin{k};
    if ischar(a)
        switch a
            case 'display'
                dodisplay = true;
            case {'singleaxes' 'multaxes'}
                axesflag = a;
            otherwise
                error('unknown flag: ''%s''',a)
        end
    elseif ishandle(a)
        ha = a;
    elseif isscalar(a)
        ncol = a;
    else
        clip = a;
    end
end
if isempty(axesflag)
    if dooutput
        axesflag = 'output';
    elseif ~isempty(ha)
        if strcmp(get(ha,'type'),'axes')
            axesflag = 'singleaxes';
        else
            axesflag = 'multaxes';
        end
    else
        axesflag = 'multaxes';
    end
end
switch axesflag
    case {'singleaxes' 'output'}
        if dodisplay
            if isempty(ha), ha = gca; end
            if ~strcmp(get(ha,'type'),'axes'), error('wrong axes handle'), end
        end
    case 'multaxes'
        hf = ha;
        if isempty(hf), hf = gcf; end
        if ~strcmp(get(hf,'type'),'figure'), error('wrong figure handle'), end
        clf(hf)
        ha = axes('parent',hf,'visible','off');
end
if isempty(clip), clip = [min(x(:)) max(x(:))]; end

% get axes size ratio
if dodisplay
    oldunit = get(ha,'units');
    set(ha,'units','pixel')
    pos = get(ha,'position');
    set(ha,'units',oldunit)
    haratio = pos(4)/pos(3);
else
    haratio = 1;
end
    
% image ratio
xratio = ny/nx;

% number of rows and columns
rowcolratio = haratio/xratio;
if any(nt==1)
    % only one condition
    nt = prod(nt);
    if isempty(ncol)
        ncol = round(sqrt(nt/rowcolratio));
        ncol = max(ncol,1);
    end
    nrow = ceil(nt/ncol);
else
    % several conditions: choose the best organization between rows and
    % coluns
    if isempty(ncol)
        testratio = nt(2)/nt(1);
        if abs(log(rowcolratio*testratio)) < abs(log(rowcolratio/testratio))
            % use nt(1) rows and nt(2) columns rather than the opposite
            nt = nt([2 1]);
            x = permute(x,[1 2 3 5 4]);
        end
    elseif ncol~=nt(1)
        if ncol==nt(2)
            nt = nt([2 1]);
            x = permute(x,[1 2 3 5 4]);
        else
            error('number of columns does not match any dimension of x')
        end
    end
    ncol = nt(1);
    nrow = nt(2);
    nt = prod(nt);
    x = reshape(x,[nx ny nc nt]);
end

% one or several axes?
if fn_ismemberstr(axesflag,{'singleaxes' 'output'})
    % make new image and fill it with frames
    if islogical(x)
        defval = false;
    else
        m = min(x(:)); M = max(x(:));
        if any(isnan(x(:)))
            defval = NaN;
        elseif m<=0 && M>=0
            defval = 0;
        elseif m<=1 && M>=1
            defval = 1;
        else
            defval = NaN;
        end
    end
    if nt<ncol*nrow, x(:,:,:,end+1:ncol*nrow) = defval; end
    x = reshape(x,nx,ny,nc,ncol,nrow);
    x = permute(x,[2 5 1 4 3]);
    if dooutput
        % add a separation between frames
        x(ny+1,:) = 0;
        x(:,:,nx+1,:) = 0;
        x = reshape(x,(ny+1)*nrow,(nx+1)*ncol,nc);
        x = [repmat(0,[1 1+(nx+1)*ncol nc]); repmat(0,[(ny+1)*nrow nc]) x]; 
    else
        x = reshape(x,ny*nrow,nx*ncol,nc);
    end
    
    % display and add lines to separate frames
    if dodisplay
        imagesc(x,'parent',ha,clip)
        axis(ha,'image')
        if ~dooutput
            for i=0:ncol
                line([1 1]*i*nx+.5,[0 ny*nrow]+.5,'color','k')
            end
            for j=0:nrow
                line([0 nx*ncol]+.5,[1 1]*j*ny+.5,'color','k')
            end
        end
    end
    
    % output?
    if dooutput
        y = permute(x,[2 1 3]);
    end
else
    % change positions of 'containing' axes to make square images
    delete(ha)
    left = pos(1); bottom = pos(2); w = pos(3); h = pos(4);
    mxratio = (ny*nrow)/(nx*ncol);
    r = haratio/mxratio; 
    if r>1
        bottom = bottom + (h-h/r)/2;
        h = h/r;
    else
        left = left + (w-w*r)/2;
        w = w*r;
    end
    w = w/ncol; h = h/nrow;
    % go! make multiple axes and display frames inside them
    k = 0;
    defunits = get(hf,'defaultaxesunits');
    for j=1:nrow
        for i=1:ncol
            k = k+1;
            if k>nt, break, end
            ha = axes('parent',hf, ...
                'units','pixel','pos',[left+(i-1)*w bottom+(nrow-j)*h w h], ...
                'units',defunits);
            imagesc(x(:,:,k)','parent',ha,clip);
            if i>1, set(ha,'yticklabel',''), end
            if j<nrow, set(ha,'xticklabel',''), end
        end
    end
end



