%% fn_flags

%% Syntax
%  [b1 ... bn] = fn_flags({val1,...,valn},{flag1,...,flagp})
%  [b1 ... bn] = fn_flags({val1,...,valn},flag1,...,flagp)
%  [b1 ... bn] = fn_flags(val1,...,valn,{flag1,...,flagp})
%  [b1 ... bn] = fn_flags(val1,...,valn,flag)
%  x = fn_flags(...)

%% Description
%  returns booleans (or vector thereof) mentioning which values appeared in
%  the list of flags, and throws an error if a flag is not in the list of
%  values 

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
