%% fn_gettrend

%% Syntax
%  [trend X beta] = fn_gettrend(y[,ind][,order][,options])

%% Description
%  Estimates a slow trend for y
%  it is possible to estimate it based on given indices only
%  
%  Inputs:
%  - y       signal
%  - ind     indices where there is supposed to be equal to baseline
%  - order   number of low-frequency regressors
%  - 'otherbase',indbis      indices where the signal is supposed to be on a
%                            plateau value
% 
%  Exemple:
% 
%  signal = [zeros(1,200) -[1:100]/100.*sin([1:100]/100*(3*pi/2)) ...
%      ones(1,300)  1+[1:100]/100.*sin([1:100]/100*(3*pi/2)) zeros(1,300)]';
%  noise = filty(1/2-rand(1000,1),200,'lm')*100;
%  y = signal + noise;
%  figure(1), plot([signal noise y])
%  legend('signal','noise','noisy signal')
%  % estimation of trend; we assume that the signal is equal to baseline for
%  % indices 1:150 and 750:1000
%  trend = fn_gettrend(y,[1:150 750:1000],6,'otherbase',350:550);
%  figure(2), plot([trend noise y-trend signal])
%  legend('estimated noise','noise','estimated signal','signal')% signal = [zeros(1,200) -[1:100]/100.*sin([1:100]/100*(3*pi/2)) ...

%% Source
% Thomas Deneux
%
% Copyright 2004-2012
%
