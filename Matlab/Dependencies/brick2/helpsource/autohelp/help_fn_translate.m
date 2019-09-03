%% fn_translate

%% Syntax
%  [y weight J dweight] = fn_translate(x,shift[,'full|valid')

%% Description
%  Translate image (or movie) x by shift.
% 
%  Input:
%  - x           2D or 3D array
%  - shift       2-element vector or 2-by-N array (N being the number of
%                frames) 
%  - shapeflag   'full' [default] or 'valid', indicate whether to return y
%                the same size as x, or only the subpart where data is
%                defined
% 
%  Output:
%  - y       the translated image or movie (points where data is not defined
%            have a NaN value)
%  - weight  image the same size of x with values btw 0 and 1: weighting
%            system for getting values of y as a continuous-derivative
%            function of shift, even at integer values of shift
%            use also logical(weight) to get the pixels whose values are
%            defined in y
%  - J       derivative of y with respect to shift (points where it is not
%            defined also have a NaN value)
%  - dweight derivative of weight with respect to shift

%% Source
% Thomas Deneux
%
% Copyright 2009-2012
%
