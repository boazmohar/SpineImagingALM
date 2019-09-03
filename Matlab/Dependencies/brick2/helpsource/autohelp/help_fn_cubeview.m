%% fn_cubeview

%% Syntax
%  [img =] fn_cubeview(data[,d[,r]])

%% Description
%  creates an image of the 3-dimensional data where we see 3 faces of a
%  cube with xy, yz and xz sides
% 
%  if no output is requested, displays the image in a figure and add edges
% 
%  faces display:       ______                                              
%      _ z            / Ypane/|                                            
%      /|           /______/Xpane                                        
%     /             |      |  |                                            
%     ---> x        |front | /                                             
%    |              |______|/                                              
%    v y                                                                   
% 
%  Input:
%  - data    3D array
%  - d       size of the side trapezes (the yz and xz faces); 
%            if d>1, size is in pixels, if d<1, it defines the ratio between
%            the side trapezes sizes and the xy size; 
%            if it is a scalar, it defines the x size, and the y size is
%            chosen so as to put perspective with and angle of 35~40
%            degrees, if it has 2 elements, it defines the x and y sizes
%  - r       ratio for perspective projection (0<r<=1, r=1 is orthogonal
%            projection)
% 
%  Output
%  - img     an image showing the cube, with a white background

%% Source
% Thomas Deneux
%
% Copyright 2010-2012
%
