function fn_saveimg(a,fname,clip,zoom,map)
% function fn_saveimg(a,fname,clip[,zoom[,cmap]])
%---
% a should be y-x-t 
% clip can be a 2-values vector, or 'fit' [default], or '?SD', or 'none'

% Thomas Deneux
% Copyright 2004-2012

if nargin<1, help fn_saveimg, return, end

if nargin<2 || isempty(fname), fname=fn_savefile; end
if nargin<3 || isempty(clip), clip='fit'; end
if nargin<4 || isempty(zoom), zoom=1; end
if nargin<5, map = []; end

% image(s) size and number color/bw
[ni nj nt nt2] = size(a);
if nt2==3
    a = permute(a,[1 2 4 3]);
    ncol = 3;
elseif nt==3
    nt = nt2;
    ncol = 3;
else 
    ncol = 1;
end

% file name
[fpath fbase fext] = fileparts(fname);
if isempty(fext)
    ext = 'png';
else
    ext = fext(2:end);
end
if nt>1
    if ~isempty(fpath), fpath = [fpath '/']; end
    fname = [fpath fbase '_'];
    lg = floor(log10(nt))+1;
    icode = ['%.' num2str(lg) 'i'];
end

% color image(s)
if ncol==3
    if ~strcmp(class(a),'uint8') && (min(a(:))<0 || max(a(:))>1)
        error('values should be between 0 and 1')
    end
    if zoom~=1
        error('no zoom allowed for color images')
    end
    a = permute(a,[2 1 3 4]); % (x,y) convention -> Matlab (y,x) convention
    if nt==1
        imwrite(a,fname,ext);
    else
        fn_progress('saving image',nt)
        for i=1:nt
            fn_progress(i)
            name = [fname num2str(i,icode) '.' ext];
            imwrite(a(:,:,:,i),name)
        end
    end
    return
end

% otherwise
a = double(a);

% clipping
if ischar(clip)
    if strcmp(clip,'fit')
        clip = [min(a(:)) max(a(:))];
    elseif findstr(clip,'SD')
        nsd = sscanf(clip,'%iSD');
        m = mean(a(:));
        sd = std(a(:));
        clip = [m-nsd*sd m+nsd*sd];
    elseif ~strcmp(clip,'none')
        error('clipping flag ''%s'' is not recognized',clip)
    end
else
    if length(clip(:))~=2
        error('clipping value is not correct')
    end
end
if ~isequal(clip,'none')
    a = (a-clip(1)) / (clip(2)-clip(1));
    a = min(max(a,0),.999);
end

% zoom parameters
if zoom~=1
    if zoom<1, disp('zoom<1 does not bin but only interpolates'), end
    if zoom>1 && mod(zoom,1)==0
        disp('integer zoom enlarges without interpolating')
        zf = true;
        ii = kron(1:ni,ones(1,zoom));
        jj = kron(1:nj,ones(1,zoom));
    else
        zf = false;
        [jj ii] = meshgrid(.5:nj-.5,.5:ni-.5);
        [jj2 ii2] = meshgrid((.5:nj*zoom-.5)/zoom,(.5:ni*zoom-.5)/zoom);
    end
end

% saving
if nt>1
    fn_progress('saving image',nt)
end
for i=1:nt
    if nt>1
        fn_progress(i)
        name = [fname num2str(i,icode) '.' ext];
    else
        name = fname;
    end
    fr = a(:,:,i)'; % (x,y) convention -> Matlab (y,x) convention
    if zoom~=1
        if zf
            fr = fr(jj,ii);
        else
            fr = interp2(jj,ii,fr,jj2,ii2,'*spline');
        end
        fr = min(max(fr,0),.999);
    end
    if ~isempty(map)
        if ischar(map), map = feval(map,256); end
        fr = floor(size(map,1)*fr)+1;
        fr = reshape(cat(3,map(fr,1),map(fr,2),map(fr,3)),nj*zoom,ni*zoom,3);
    else % gray image
        %fr = floor(length(map)*fr)+1;
    end
    imwrite(fr,name,ext)
end

