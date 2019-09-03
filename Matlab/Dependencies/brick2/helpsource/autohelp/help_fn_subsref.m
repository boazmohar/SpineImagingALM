%% fn_subsref

%% Syntax
%  B = fn_subsref(A,idx1,idx2,...)
%  indices = fn_subsref(siz,idx1,idx2,...[,'global|local'])

%% Description
%  
%  Input:
%  - A/siz       array, or size of an array (the function guesses that it is
%                a size vector if the argument is a vector of size less than
%                or equal to 5, therefore, an array A should be larger than
%                this)
%  - idx1, ...   indices in the successive dimensions; can be either numeric
%                values (the indices themselves) or strings that are
%                interpreted (such as ':', '1:3', '1 3:4', etc.)
%                the number of indices specification must be one or the
%                number of dimensions of A / the length of siz
%  - 'global|local'  indicates whether the output indices should be a vector
%                of global indices [default behavior], or a cell array of
%                indices for each coordinates
% 
%  Output:
%  - B/indices   sub-array formed by the elements of A specified by the
%                indices idx1, ..., or indices formatted as 'global' or
%                'local'
% 
%  Examples:
%  - fn_subsref([1 2 3; 4 5 6],':',2:3)   returns the sub-array [2 3; 5 6]
%  - fn_subsref([1 2 3; 4 5 6],'1','1 3') returns the sub-array [1 3]
%  - fn_subsref([2 3],'1','1 3','local')  returns the indices {[1] [1 3]}
%  - fn_subsref([2 3],'1','1 3')          returns the indices [1 5]

%% Source
% Thomas Deneux
%
% Copyright 2006-2012
%
