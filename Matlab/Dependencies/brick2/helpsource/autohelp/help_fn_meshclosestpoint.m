%% fn_meshclosestpoint

%% Syntax
%  [i p] = fn_meshclosestpoint(vert,p)
%  [i p] = fn_meshclosestpoint(mesh)

%% Description
%  Input:
%  - vert    (3xN) array of points 
%            (can be alternatively a {vertices,faces} array, then first cell
%            array element is used)
%  - p       3D point (vector) or 3D line (2x3 array as the output of
%            get(gca,'CurrentPoint'))
%  - mesh    {vertices,faces} : mesh is displayed and then function waits for
%            button press
% 
%  Output:
%  - i       indice of the vertex that is closest to p

%% Source
% Thomas Deneux
%
% Copyright 2005-2012
%
