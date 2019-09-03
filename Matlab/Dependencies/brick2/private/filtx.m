function y = filtx(x, s, type, s1)
%FILTX	Gaussian filter for rows of a matrix.
%	Y = FILTX(X, S, TYPE, S1) 
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
%	See also: FILTY, FILT2, MFILTX, MFILTY, MFILT2
%
%	Doron, 21 Aug 1994.
%	Last revision: 25 Apr 1996.

if (nargin < 2), s = 1; end;
if (nargin < 3), type = 'l'; end;
if (nargin < 4), s1 = 1.5 * s; end;

mind = find(type == 'm');
do_mirror = length(mind);
if do_mirror,
	type(mind) = [];
	[m,n,dum]=size(x);
    % Padding with 4*sigma pixels.
    if ismember('b',type), pad = ceil(4*s1); else pad = ceil(4*s); end
	pn = min(pad, n);	% Pad is not more than image width
	x = x(:,[pn:-1:1 1:n n:-1:(n-pn+1)],:);
	y = filtx(x, s, type, s1);
	y = y(:,(1:n)+pn,:);	% get back to original size
	return;
end;

doz = any(type=='z');
type(type=='z')=[];

mn=size(x);
sf  = 1 / (2*pi*s);	% sigma in frequency space
sf1 = 1 / (2*pi*s1);	%

%xx = fftshift([-n/2:1:(n/2-1)]')/n;
% The above is correct only for even n. The next line is more general.
xx = [0:ceil(mn(2)/2-1), ceil(-mn(2)/2):-1]' / mn(2);

g = exp(-xx.^2/(2*sf.^2));
switch type
    case 'h'
        g = 1-g;
    case 'b'
        g1 = exp(-xx.^2/(2*sf1.^2));
        g = g - g1;
    case 's'
        g1 = exp(-xx.^2/(2*sf1.^2));
        g = 1 - (g-g1);
end
if doz
    g(1) = 1-g(1);
end
g = repmat(g(:),[1 mn(1) mn(3:end)]);


fx = fft(permute(x,[2 1 3:ndims(x)]));
y = permute(real(ifft(fx.*g)),[2 1 3:ndims(x)]);


