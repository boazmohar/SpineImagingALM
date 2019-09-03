%% fn_meshselectpoint

%% Syntax
%  [i p tr] = fn_meshselectpoint(mesh[,p])

%% Description
%  Input:
%  - mesh    {vertices,faces}
%  - p       3D line (2x3 array as the output of get(gca,'CurrentPoint'))
%  - mesh    {vertices,faces} : mesh is displayed and then function waits for
%            button press
% 
%  Output: 
%  - i       indice of the selected vertex
%  - p       coordinate of that vertex
%  - tr      indice of the selected triangle

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
