function a = fn_readmovie(filename)
% function a = fn_readmovie(filename)
%---
% read an avi file and stores it into a 2D+time array 
% (2D+time+channel if color movie)
%
% See also fn_savemovie

% Thomas Deneux
% Copyright 2004-2012

if nargin<1
    filename = fn_getfile('*.avi');
end

a = aviread(filename);

switch size(a(1).cdata,3)
    case 1
        a = cat(3,a.cdata);
    case 3
        s = size(a(1).cdata);
        nt = length(a);
        a = cat(2,a.cdata);
        a = reshape(a,[s(1) s(2) nt 3]);
    otherwise
        error('problem')
end 