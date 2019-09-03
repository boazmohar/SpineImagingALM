%% fn_framedesign

%% Syntax
%  [pos cmd] = fn_framedesign([grob[,pos|cmd[,resetflag]]])

%% Description
%  utility for positioning frames in a figure
% 
%  Input:
%  - grob        structure with fields the names of graphic objects
%                considered and values their handles; their must be one
%                field 'hf' for the main figure
%  - pos         structure with fields the names of graphic objects and
%                values their positions; fields are allowed to differ from
%                those of 'grob'
%  - cmd         string defining pos (should be something like
%                'pos.hf=[..]; ...')
%  - resetflag   true [default] or false - user resize?
% 
%  Output:
%  - pos         as above
%  - cmd         as above

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
