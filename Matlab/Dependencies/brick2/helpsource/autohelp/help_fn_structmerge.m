%% fn_structmerge

%% Syntax
%  s = fn_structmerge(s,s1[,'skip|strict'][,'recursive'][,'type'])

%% Description
%  set or replace values in s from those in s1, where s and s1 are
%  structures of the same size
%  - if 'skip', or 'strict' flag is specified: does not add new field in
%    structure s (generates error if 'strict' flag and s1 has additional
%    fields)
%  - if 'recursive' flag, field values which are themselves structures are
%    not merely replaced, but are also merged using fn_structmerge 
%  - if 'type' flag is specified: also requires field values to be the same
%    class in s and s1 when the field already exists in s (generates error
%    if it is not the case, except that it performs the conversions
%    0/1->false/true  and char array->cell array of strings)

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
