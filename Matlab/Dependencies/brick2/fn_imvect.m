function x = fn_imvect(x,mask,varargin)
% function x = fn_imvect(x,mask[,outputtype][,outsidevalue])
% function k = fn_imvect(ij,mask[,outsidevalue])
% function ij = fn_imvect(k,mask[,outsidevalue])
%---
% switch between "image" and "vector" representation of the pixels in an
% image
%
% Input:
% - x       array of size (nx,ny,nt,...) or (np,nt,...)
% - mask    logical array of size (nx,ny), such that sum(mask(:))==np
% - outputtype      'vector' or 'image': default behavior toggles
%                   represenation, by setting outputtype, x is unchanged if
%                   is already has the desired representation
% - outsidevalue    value to set outside the mask in the image [default=0]
% 
% Output:
% - x       array size became (np,nt,...) or (nx,ny,nt,...), respectively
%
% See also fn_indices

% Thomas Deneux
% Copyright 20011-2012

% input
if ~islogical(mask), error('mask must be a logical array'), end
outsidevalue = 0;
outputtype = 'toggle';
for k=1:nargin-2
    a = varargin{k};
    if isnumeric(a)
        outsidevalue = a;
    else
        if ~fn_ismemberstr(a,{'image' 'vector' 'toggle'}), error('invalid flag ''%s''',a), end
        outputtype = a;
    end
end

[nx ny] = size(mask);
np = sum(mask(:));
if isvector(x), x = x(:); end
s = size(x);
if all(s(1:2)==[nx ny])                     % image to vector
    if strcmp(outputtype,'image'), return, end % (nothing to do)
    s1 = s(3:end);
    x = reshape(x,[nx*ny s1 1]);
    x = x(mask,:);
    x = reshape(x,[np s1 1]);
elseif s(1)==np                             % vector to image
    if strcmp(outputtype,'vector'), return, end % (nothing to do)
    s1 = s(2:end);
    xold = x;
    if outsidevalue==0
        x = zeros([nx*ny s1],class(x));
    elseif outsidevalue==1
        x = ones([nx*ny s1],class(x));
    elseif isnan(outsidevalue)
        x = nan([nx*ny s1],class(x));
    else
        x = outsidevalue*ones([nx*ny s1],class(x));
    end
    x(mask,:) = xold(:,:);
    x = reshape(x,[nx ny s1]);
elseif isscalar(x)                          % vector index to image index
    k = x;
    test = zeros(np,1); test(k) = 1;
    test = fn_imvect(test,mask);
    [i j] = find(test);
    x = [i j];
elseif isvector(x) && length(x)==2          % image index to vector index
    ij = x;
    i = ij(1); j = ij(2);
    test = false(nx,ny); test(i,j) = true;
    test = fn_imvect(test,mask);
    k = find(test);
    x = k;
else
    error('dimensions of x and mask do not fit')
end
