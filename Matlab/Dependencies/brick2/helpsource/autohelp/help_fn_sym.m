%% fn_sym

%% Syntax
%  s = fn_sym(a[,uplo])
%  a = fn_sym(s[,uplo])
%  idx = fn_sym(ij[,uplo])
%  idx = fn_sym(i,j[,uplo])
%  [i j] = fn_sym(idx,[,uplo])
%  ij = fn_sym(idx[,uplo])

%% Description
%  converts square (symmetric) matrix to a compact vector
%  or inverse conversion
% 
%  uplo : 'U' [default] or 'L'
%  use 'U' to match the C++ implementation

%% Source
% Thomas Deneux
%
% Copyright 2003-2012
%
