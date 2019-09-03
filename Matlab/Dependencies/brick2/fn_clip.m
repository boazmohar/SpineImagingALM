function x = fn_clip(x,clipflag,outflag)
% function x = fn_clip(x,clipflag,outflag)
%---
% Rescale and restrict to a specific range ("clip") the data in an array.
% Make a color image if requested.  
%
% Input:
% - x           array (any dimension)
% - clipflag    clipping mode:
%               [a b]                   define manually min and max value
%               'fit','mM' or 'minmax'  use minimum and maximum [default]
%               'Xstd'                  use mean and X times standard deviation
% - outflag     output format
%               [a b]       define minimum and maximum value [default, with
%                           a=0 and b=1]
%               n           integer values between 1 and n
%               nx3 array   returns a (n+1)-dimensional array using this
%                           colormap 
%               char array  use this colormap (for example 'jet' -> use
%                           jet(256))

% Thomas Deneux
% Copyright 2007-2012

if nargin==0, help fn_clip, return, end

x = double(x);
if nargin<2 || isempty(clipflag), clipflag='mM'; end
if nargin<3, outflag=[0 1]; end

% clipping mode
if isnumeric(clipflag)
    if ~isvector(clipflag) || length(clipflag)~=2, error('clipping vector must have 2 elements'), end
    clip = clipflag;
else
    xstd = regexpi(clipflag,'^([\d.]*)st{0,1}d$','tokens');
    if ~isempty(xstd)
        xstd = xstd{1}{1};
        if isempty(xstd), xstd=1; else xstd=str2double(xstd); end
        m = mean(x(:));
        st = std(x(:));
        clip = m + [-1 1]*xstd*st;
    elseif fn_ismemberstr(clipflag,{'fit' 'mM' 'minmax'})
        clip = [min(x(:)) max(x(:))];
    else
        error('erroneous clipping option')
    end
end

% clip
x = (x-clip(1))/diff(clip);
x = min(1-eps(single(1)),max(0,x)); % it is convenient that 1 cannot be reached

% output mode
if ischar(outflag)
    docolor = true;
    cm = feval(outflag,256);
    n = size(cm,1);
elseif ~isvector(outflag)
    if size(outflag,2)~=3, error('colormap must have 3 columns'), end
    docolor = true;
    cm = outflag;
    n = size(cm,1);
elseif isscalar(outflag)
    docolor = false;
    n = outflag;
    if mod(n,1) || n<=0, error('scalar for output format must be a positive integer'), end
else
    docolor = false;
    if length(outflag)~=2, error('vector for output format must have 2 elements'), end
    n = 0;
    a = outflag(1);
    b = outflag(2);
end

% scaling
if n
    x = floor(n*x)+1; % n+1 cannot be reached    
else
    x = a + x*(b-a); % b cannot be reached
end

% color
if docolor
    s = size(x); if s(2)==1, s(2)=[]; end
    x = reshape(cm(x,:),[s 3]);
    if length(s)>2, x = permute(x,[1 2 length(s)+1 3:length(s)]); end
end




        