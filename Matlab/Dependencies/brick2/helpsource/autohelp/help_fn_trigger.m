%% fn_trigger

%% Syntax
%  y = fn_trigger(x,indices[,dim][,nframes])

%% Description
%  average according to triggers
% 
%  Input:
%  - x           data to trigger
%  - indices     trigger position, in indices coordinates
%  - dim         dimension of the data on which to trigger and average [by
%                default, the last dimension is used]
%  - nframes     size of the output in this dimension: it should be a
%                2-element vector (how many frames before trigger, and how
%                many frames after trigger), or a scalar (the same number is
%                used both times) [the default is 0, i.e. only averaging at
%                the trigger positions]

%% Source
% Thomas Deneux
%
% Copyright 2008-2012
%
