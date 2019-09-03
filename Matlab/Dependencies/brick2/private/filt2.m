function y = filt2(x, s, type, s1)
%FILT2	2D Gaussian filter
%	Y = FILT2(X, S, TYPE, S1) 
%
%	X:	Matrix to be filtered
%	S:	Sigma of gaussian to use (default 1)
%	TYPE:	A string of 1-2 characters to determind type of filter.
%		'l' (default) for low-pass, 'h' for high-pass,
%		'b' for band-pass (difference of gaussians)
%		If also 'm' is given (e.g, 'lm'), the image will be padded
%		with it's mirror reflections before filtering, instead of the
%		default which effectively does wrap around.
%	S1:	Sigma for high pass in case of band-pass (default 1.5*S)
%
%	See also: FILTX, FILTY, MFILT2, MFILTX, MFILTY
%
%	Doron, 21 Aug 1994.
%	Last revision: 25 Apr 1996.

if (nargin < 2) s = 1; end;
if (nargin < 3) type = 'l'; end;
if (nargin < 4), s1 = 1.5 * s; end;


% start - added by T. Deneux
if ndims(x)>=3
    y = x;
    [nx ny n3] = size(x);
    for k=1:n3, y(:,:,k) = filt2(x(:,:,k),s,type,s1); end
    return
end
% end   - added by T. Deneux
        
mind = find(type == 'm');
do_mirror = length(mind);
if do_mirror,
	type(mind) = [];
	[m,n]=size(x);
    % Padding with 4*sigma pixels.
    if ismember('b',type), pad = ceil(4*s1); else pad = ceil(4*s); end
	pm = min(pad, m);	% Pad is not more than image size
	pn = min(pad, n);
	x = x([pm:-1:1 1:m m:-1:(m-pm+1)], [pn:-1:1 1:n n:-1:(n-pn+1)]);
	y = filt2(x, s, type, s1);
	y = y([1:m]+pm, [1:n]+pn);	% get back to original size
	return;
end;

[m,n]=size(x);
sf  = 1 / (2*pi*s);	% sigma in frequency space
sf1 = 1 / (2*pi*s1);	%

%[xx,yy]= meshgrid(fftshift([-n/2:1:(n/2-1)]/n), fftshift([-m/2:1:(m/2-1)]/m));
% The above is correct only for even m and n. The next 3 lines are more general
xind = [0:ceil(n/2-1), ceil(-n/2):-1] / n;
yind = [0:ceil(m/2-1), ceil(-m/2):-1] / m;
[xx,yy]= meshgrid(xind, yind);

g = exp(-(xx.^2+yy.^2)/(2*sf.^2));
if ismember('z',type), g(1,1) = 0; end
if ismember('h',type), g = 1-g; end;
if ismember('b',type),
	g1 = exp(-(xx.^2+yy.^2)/(2*sf1.^2));
	g = g - g1;
    if ismember('z',type), g(1,1) = 1; end
end;

fx = fft2(x);

y = real(ifft2(fx.*g));
