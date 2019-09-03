%% fn_timevector

%% Syntax
%  times|count = fn_timevector(times|count,dt|tidx[,outputtype])

%% Description
%  switch representation of a point process between the times of the events
%  and the number of events per time bin
% 
%  Input:
%  - count       vector of spike count at each time bin (or array or cell
%                array of vectors)
%  - times       vector of spike times (or cell array thereof)
%  - dt          scalar - size of time bin
%  - tidx        vector of all time instants
%  - outputtype  'times', 'rate' or 'count' [default: output type is the
%                opposite of input type]; use 'rateperperiod' or
%                'countperperiod' to count using time bins that are not
%                centered on the time points in tidx, but rather whose sides
%                are defined by the time points in tidx (in such case, the
%                number of bins will be one less than the number of time
%                points)
% 
%  Output:
%  - count       column vector or array
%  - times       row vector or cell array of row vectors

%% Source
% Thomas Deneux
%
% Copyright 2010-2012
%
