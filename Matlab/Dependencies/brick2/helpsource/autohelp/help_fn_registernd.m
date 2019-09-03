%% fn_registernd

%% Syntax
%  [mov check] = fn_registernd(ref,test,sigma,hp)

%% Description
%  register z-stack test to the z-stack ref and returns the distance to move
%  in pixels
%  
%  Input:
%  - ref         ND array
%  - test        ND array, sizes smaller than ref
%  - sigma       nrep*ndim array: value of low-pass filter to apply first
%                [default: 1]
%  - hp          parent for display: either figure or uipanel handle
%                [default: no display]
%  
%  Output:
%  - mov         estimated distance in pixels to go from the center of ref
%                to the place where the center of test has been register
%  - check       interpolation of ref at the registered location, the same
%                size of test, allows to check that the alignment is correct
% 
%  the maximal distance is such that the portion of 'test' which is allowed to
%  move out of 'ref' is less than 1/4 in each dimension

%% Source
% Thomas Deneux
%
% Copyright 20011-2012
%
