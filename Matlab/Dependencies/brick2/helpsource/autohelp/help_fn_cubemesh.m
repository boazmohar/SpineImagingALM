%% fn_cubemesh

%% Syntax
%  fn_cube(data[,fname])

%% Description
%  creates data mesh and its texture showing the 6 faces of data 3D-object
%  if called without ouptput arguments, displays it in data figure, or save it
%  in data file 
%  
%  faces display:       ______                                              
%      _ z            /  E   /| C                                          
%      /|        D  /______/  |                                           
%     /             |      | B|                                            
%     ---> x        | A    | /                                             
%    |              |______|/                                              
%    v y                F                                                  
% 
%  Input:
%  - data       3D array
%  - fname   file name with no extension: automatically, extensions .tri and
%            .text will be added for the mesh and texture files
% 
%  Output: 
%  - m       cell array {vertices faces} representing the mesh
%  - text    vector with values for each vertex

%% Source
% Thomas Deneux
%
% Copyright 2010-2012
%
