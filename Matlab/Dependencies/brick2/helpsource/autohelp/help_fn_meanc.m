%% fn_meanc

%% Syntax
%  [m, lb, ub] = fn_meanc(x[,dim[,alpha]])

%% Description
%  same as mean, but also returns a alpha level confidence interval for the mean
%  default dim is 1, default alpha is 90%
%  si mu vraie moyenne, m moyenne estim�e, s variance estim�e (sans biais)
%  alors (m-mu)*sqrt(n)/s ~ Student(n-1)

%% Source
% Thomas Deneux
%
% Copyright 2004-2012
%
