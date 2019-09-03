%% fn_linespecs

%% Syntax
%  lineoptions = fn_linespecs(line options as in plot)
%  lineoptions = fn_linespecs({line options as in plot})

%% Description
%  If necessary, converts the first argument into the accurate set of line
%  properties {Propertyname1,Value1,...} according to the same syntax as
%  the plot command.
%  This can be used to set line properties more easily.
%  ex: opt=fn_linespecs('r:','linewidth',2);
%      line([0 1],[0 1],opt{:})
%  
%           b     blue          .     point              -     solid
%           g     green         o     circle             :     dotted
%           r     red           x     x-mark             -.    dashdot 
%           c     cyan          +     plus               --    dashed   
%           m     magenta       *     star             (none)  no line
%           y     yellow        s     square
%           k     black         d     diamond
%           w     white         v     triangle (down)
%                               ^     triangle (up)
%                               <     triangle (left)
%                               >     triangle (right)
%                               p     pentagram
%                               h     hexagram

%% Source
% Thomas Deneux
%
% Copyright 2004-2012
%
