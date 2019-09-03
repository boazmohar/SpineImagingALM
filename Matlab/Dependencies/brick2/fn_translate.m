function [y weight J dweight] = fn_translate(x,shift,shapeflag)
% function [y weight J dweight] = fn_translate(x,shift[,'full|valid')
%---
% Translate image (or movie) x by shift.
%
% Input:
% - x           2D or 3D array
% - shift       2-element vector or 2-by-N array (N being the number of
%               frames) 
% - shapeflag   'full' [default] or 'valid', indicate whether to return y
%               the same size as x, or only the subpart where data is
%               defined
%
% Output:
% - y       the translated image or movie (points where data is not defined
%           have a NaN value)
% - weight  image the same size of x with values btw 0 and 1: weighting
%           system for getting values of y as a continuous-derivative
%           function of shift, even at integer values of shift
%           use also logical(weight) to get the pixels whose values are
%           defined in y
% - J       derivative of y with respect to shift (points where it is not
%           defined also have a NaN value)
% - dweight derivative of weight with respect to shift

% Thomas Deneux
% Copyright 2009-2012

% Input
if nargin<3, shapeflag='full'; end

% Bi-cubic interpolation needs a grid of 4x4 points to define value in one
% point, i.e. 2 values in each left/right/up/down direcion.
% If we want to shift the image rightward, we need to look leftward when
% filling the new image 
dec = -double(shift);
% first attempt is deci = floor(dec), but for example if dec = 5, it is
% possible and preferrable to have deci = 4, because this will allows
% getting one more value defined on the right, hence the complicated
% formula below
deci = zeros(1,2);
for k=1:2
    if dec(k)>0, deci(k)=ceil(dec(k))-1; else deci(k)=floor(dec(k)); end
end
% coordinates inside the grid
u = dec(1)-deci(1);
v = dec(2)-deci(2);

% Parts of the input and output images to consider
[ni nj nt] = size(x);
iin = max(1+(deci(1)-1),1):min(ni+(deci(1)+2),ni);
jin = max(1+(deci(2)-1),1):min(nj+(deci(2)+2),nj);
iout = max(1,1-(deci(1)-1)):min(ni,ni-(deci(1)+2));
jout = max(1,1-(deci(2)-1)):min(nj,nj-(deci(2)+2));
% corresponding weight
if nargout>=2
    u2 = max(eps,min(1-eps,3*u^2-2*u.^3));
    v2 = max(eps,min(1-eps,3*v^2-2*v^3));
    weight = zeros(ni,nj);
    weight(iout,jout) = 1;
    dod = (nargout>=3);
    if dod, dweight = zeros(ni,nj,2); end
    if iout(1)  == 1-(deci(1)-1)
        if dod
            dweight(iout(1),:,1)   = weight(iout(1),:)*(-6*u*(1-u)); 
        end
        weight(iout(1),:)   = weight(iout(1),:)*u2; 
    end
    if iout(end)== ni-(deci(1)+2)
        if dod
            dweight(iout(end),:,1)  = weight(iout(end),:)*(+6*u*(1-u));
        end
        weight(iout(end),:) = weight(iout(end),:)*(1-u2);     
    end
    if jout(1)  == 1-(deci(2)-1)
        if dod
            dweight(:,jout(1),1)   = dweight(:,jout(1),1)*v2; 
            dweight(:,jout(1),2)   = weight(:,jout(1))*(-6*v*(1-v)); 
        end 
        weight(:,jout(1))   = weight(:,jout(1))*v2; 
    end
    if jout(end)== nj-(deci(2)+2)
        if dod
            dweight(:,jout(end),1) = dweight(:,jout(end),1)*(1-v2);
            dweight(:,jout(end),2) = weight(:,jout(end))*(+6*v*(1-v));
        end 
        weight(:,jout(end)) = weight(:,jout(end))*(1-v2); 
    end
    total  = sum(weight(:));
    if dod
        dtotal = sum(sum(dweight));
        dweight = dweight/total - cat(3,weight*dtotal(1),weight*dtotal(2))/total^2;
    end
    weight  = weight/total;
end
% take the needed part of x
y0 = x(iin,jin);

% Spline Interpolation!
% First the matrix from http://en.wikipedia.org/wiki/Bicubic_interpolation
A = 1/2 * [0 2 0 0; -1 0 1 0; 2 -5 4 -1; -1 3 -3 1];
% Convolution in the first dimension
u = dec(1)-deci(1);
filti = [1 u u^2 u^3]*A;
filti = flipud(filti'); % using convolution instead of correlation: must flip the filter
y1 = convn(y0,filti,'valid');
% Convolution in the second dimension
v = dec(2)-deci(2);
filtj = [1 v v^2 v^3]*A;
filtj = fliplr(filtj); 
y12 = convn(y1,filtj,'valid');

% Derivatives!
if nargout>=3
    if nt>1, error('cannot compute derivative for a movie'), end
    dfiltidu = [0 1 2*u 3*u^2]*A;
    dfiltidu = flipud(dfiltidu');
    dfiltidv = [0 1 2*v 3*v^2]*A;
    dfiltidv = fliplr(dfiltidv);
    y2 = convn(y0,filtj,'valid');
    dy12du = convn(y2,dfiltidu,'valid');
    dy12dv = convn(y1,dfiltidv,'valid');
    dy12dshift = cat(3,-dy12du,-dy12dv);
end
    
% Fill-in the output
switch shapeflag
    case 'full'
        y = nan(ni,nj,nt);
        y(iout,jout,:) = y12;
        if nargout>=3
            J = nan(ni,nj,2);
            J(iout,jout,:) = dy12dshift;
        end
    case 'valid'
        y = y12;
        if nargout>=3, J = dy12dshift; end
    otherwise
        error('unknown shape flag ''%s''',shapeflag)
end


