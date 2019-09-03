function y = fn_meantc(x)
% function y = fn_meantc(x)
%---
% returns the average time course of a movie

% Thomas Deneux
% Copyright 2006-2012

s = size(x);
x = reshape(x,[s(1)*s(2) s(3:end)]);
y = mean(x,1);
y = shiftdim(y,1);