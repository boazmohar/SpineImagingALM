%% fn_clip

%% Syntax
%  x = fn_clip(x,clipflag,outflag)

%% Description
%  Rescale and restrict to a specific range ("clip") the data in an array.
%  Make a color image if requested.  
% 
%  Input:
%  - x           array (any dimension)
%  - clipflag    clipping mode:
%                [a b]                   define manually min and max value
%                'fit','mM' or 'minmax'  use minimum and maximum [default]
%                'Xstd'                  use mean and X times standard deviation
%  - outflag     output format
%                [a b]       define minimum and maximum value [default, with
%                            a=0 and b=1]
%                n           integer values between 1 and n
%                nx3 array   returns a (n+1)-dimensional array using this
%                            colormap 
%                char array  use this colormap (for example 'jet' -> use
%                            jet(256))

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
