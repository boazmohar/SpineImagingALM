%% fn_eegplot

%% Syntax
%  hl = fn_eegplot(usual plot arguments[,stepflag]['num'])
% Like Matlab  'plot', but separates line by a distance 'ystep'

%% Description
%  the way 'ystep' is calculated is determinated by 'stepflag':
%    - x         use ystep = x
%    - 'STD'     use ystep = mean(std(data))
%    - 'xSTD'    use ystep = x*mean(std(data)) [default is STD]
%    - 'fit' 	use ystep = max(max(data)-min(data))
%    - 'xfit'    use ystep = x * max(max(data)-min(data))
% 
%  additional flag 'num' indicates to rescale the data so that numbers on
%  the ordinates axis correspond to data number

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
