function varargout = fn_register(varargin)
% function [shift e xreg] = fn_register(x,par)
% function par = fn_register('par')

% Thomas Deneux
% Copyright 2011-2012

if nargin==0, help fn_register, end

x = varargin{1};
if ischar(x)
    if ~strcmp(x,'par'), error argument, end
    varargout = {defaultpar};
else
    par = defaultpar;
    if nargin>=2, par = fn_structmerge(par,varargin{2},'skip'); end
    nout = max(nargout,1);
    varargout = cell(1,nout);
    [varargout{:}] = register(x,par);
end
%if nargout==0, varargout = {}; end

%---
function par = defaultpar

par.maxshift = .5;      % shift cannot exceed half of the image size [TODO!]
par.ref = 10;           % reference frame: use average of first 10 frames as default
par.repeat = false;     % repeat the estimation with the average resampled movie as new reference frame
par.smoothing = 0;      % smoothing of the estimated drift!!
par.display = 'none';   % possibilities: 'iter' and 'final'
par.maxerr = [];        % maximal error allowed: if not attained, try other initializations; [TODO!]
par.shift0 = [0 0];

%---
function [shift e xreg] = register(x,par)

% Size and indices of reference image that will always stay inside the
% image
x = double(x);
[ni nj nt] = size(x);
[par.ni par.nj par.nt] = size(x);

% Reference frame
if isscalar(par.ref)
    ref = mean(x(:,:,1:min(par.ref,nt)),3);
else
    ref = par.ref;
    if any(size(ref)~=[ni nj]), error 'size mismatch with reference frame', end
end
ref = (ref-mean(ref(:)))/std(ref(:)); % normalize image

% Maximal move
if isscalar(par.maxshift)
    if par.maxshift>1
        par.maxshift = par.maxshift*[1 1];
    else
        par.maxshift = par.maxshift*[ni nj];
    end
end

% Display
if ~strcmp(par.display,'none')
    par.hf = figure(678);
    set(par.hf,'numberTitle','off','name','fn_register')
end

% Register
opt = optimset('Algorithm','active-set','GradObj','off', ...
    'tolx',1e-4,'tolfun',1e-4,'maxfunevals',1000,'display','none');
shift = zeros(2,nt);
if nt>1, fn_progress('register frame',nt), end
d = par.shift0(:);
if nargout>=2, e = zeros(1,nt); end
for k=1:nt
    if nt>1 && ~mod(nt-k,10), fn_progress(k), end
    xk = x(:,:,k);
    xk = (xk-mean(xk(:)))/std(xk(:)); % normalize image
    d = fmincon(@(d)energy(d,xk,ref,par),d, ...
        [],[],[],[],-par.maxshift,par.maxshift,[],opt);
    if strcmp(par.display,'final')
        figure(par.hf), fn_alignimage(x',ref',d([2 1]))
    end
    shift(:,k) = d;
    if nargout>=2, e(k) = energy(d,xk,ref,par); end
end

% Smooth estimated drift
if par.smoothing
    shift = [repmat(shift(:,1),1,nt) shift repmat(shift(:,nt),1,nt)];
    shift = filtx(shift,par.smoothing);
    shift = shift(:,nt+1:2*nt);
end
if nargout<3, return, end

% Resample
if nt>1, disp('resample frames'), end
xreg = zeros(ni,nj,nt);
for k=1:nt
    xreg(:,:,k) = fn_translate(x(:,:,k),-shift(:,k),'full');
end



%---
function [e de] = energy(d,x,ref,par)

doJ = (nargout==2);
if ~doJ
    [xpred weight] = fn_translate(ref,d,'valid');
else
    [xpred weight J dweight] = fn_translate(ref,d,'valid');
end
mask = logical(weight);
xpred = xpred(:);
N = length(xpred);
dif = xpred - x(mask);
dif2 = dif.^2;
weight = weight(mask);

e  = sum(dif2.*weight);
if doJ
    % TODO: still, something is not good with the derivative, in particular at integer values of shift
    J = reshape(J,[N 2]); 
    dweight = reshape(dweight,[par.ni*par.nj 2]);
    dweight = dweight(mask,:);
    
    de = 2*(dif.*weight)'*J + sum(repmat(dif2,1,2).*dweight);
end

if strcmp(par.display,'iter')
    figure(par.hf), fn_alignimage(x',ref',d([2 1]))
end

