function y = filty(x, s, type, s1)
%FILTY	Gaussian filter for columns of a matrix.
%	Y = FILTY(X, S, TYPE, S1)
%
%	X:	Matrix to be filtered
%	S:	Sigma of gaussian to use (default 1)
%	TYPE:	A string of 1-2 characters to determind type of filter.
%		'l' (default) for low-pass, 'h' for high-pass,
%		'b' for band-pass (difference of gaussians)
%       's' for stop-pass
%		If also 'm' is given (e.g, 'lm'), the image will be padded
%		with it's mirror reflections before filtering, instead of the
%		default which effectively does wrap around.
%       [Thomas Deneux:] If also 'z' is given (like 'zero'), it also
%       removes the average (in the case of low-pass), or on the contrary
%       it does not remove the average (in the case of high-pass and
%       band-pass)
%	S1:	Sigma for high pass in case of band-pass (default 1.5*S)
%
%	See also: FILTX, FILT2, MFILTX, MFILTY, MFILT2
%
%	Doron, 21 Aug 1994.
%	Last revision: 25 Apr 1996. + T. Deneux, later

if (nargin < 2) s = 1; end;
if (nargin < 3) type = 'l'; end;
if (nargin < 4), s1 = 1.5 * s; end;

if ~ischar(type) || ~isvector(type)
	error('TYPE should be a 1- to 3-element string vector!');
end;

siz = size(x);
[m,n]=size(x);
do_mirror = ismember('m',type);
do_imirror = ismember('i',type);
do_circular = ismember('c',type);
do_replicate = ismember('r',type);
type = intersect(type,'lhbsz');
if do_mirror || do_imirror || do_replicate,
    % Padding with 4*sigma pixels.
    if ismember('b',type), pad = ceil(4*s1); else pad = ceil(4*s); end
	pm = min(pad, m);	% Pad is not more than image height
    if do_mirror
    	x = x([pm:-1:1 1:m m:-1:(m-pm+1)], :);
    elseif do_imirror % inverse mirror -> will force y(1)=0 and y(end)=x(end)
        x = [-x(pm:-1:1,:) ; x(1:m,:) ; (2*x(m)-x(m:-1:(m-pm+1),:))];
    elseif do_circular
        x = x([m-pm+1:m 1:m 1:pm], :);
    else % replicate extremities
        x = [x(ones(pm,1),:) ; x(1:m,:) ; x(m*ones(pm,1),:)];
    end
	y = filty(x, s, type, s1);
	y = y((1:m)+pm, :);	% get back to original size
    y = reshape(y,siz);
	return;
end;
do_zero = ismember('z',type);
type = intersect(type,'lhbs');

sf  = 1 / (2*pi*s);	% sigma in frequency space
sf1 = 1 / (2*pi*s1);	%

%yy = fftshift(ceil([-m/2:1:(m/2-1)]'))/m;
% The above is correct only for even m. The next line is more general.
yy = [0:ceil(m/2-1), ceil(-m/2):-1]' / m;

switch type
    case 'l'
        g = exp(-yy.^2/(2*sf.^2));
    case 'h'
        g = 1-exp(-yy.^2/(2*sf.^2));
    case 'b'
        g = exp(-yy.^2/(2*sf.^2)) - exp(-yy.^2/(2*sf1.^2));
    case 's'
        g = 1-exp(-yy.^2/(2*sf.^2)) + exp(-yy.^2/(2*sf1.^2));
end;

if do_zero
    g(1) = 1-g(1);
end

g = g(:,ones(1,n));
g = reshape(g,size(x));

fx = fft(x);

y = real(ifft(fx.*g));
y = reshape(y,siz);
