%% fn_imvect

%% Syntax
%  x = fn_imvect(x,mask[,outputtype][,outsidevalue])
%  k = fn_imvect(ij,mask[,outsidevalue])
%  ij = fn_imvect(k,mask[,outsidevalue])

%% Description
%  switch between "image" and "vector" representation of the pixels in an
%  image
% 
%  Input:
%  - x       array of size (nx,ny,nt,...) or (np,nt,...)
%  - mask    logical array of size (nx,ny), such that sum(mask(:))==np
%  - outputtype      'vector' or 'image': default behavior toggles
%                    represenation, by setting outputtype, x is unchanged if
%                    is already has the desired representation
%  - outsidevalue    value to set outside the mask in the image [default=0]
%  
%  Output:
%  - x       array size became (np,nt,...) or (nx,ny,nt,...), respectively
% 
%  See also fn_indices

%% Source
% Thomas Deneux
%
% Copyright 20011-2012
%
