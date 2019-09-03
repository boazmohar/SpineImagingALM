function a = fn_readimg(fname,flag)
% function a = fn_readimg(fname[,'permute'])
%---
% read image using imread, and handles additional features:
% - converts to double
% - detects if color or gray-scale images (in the last case, use a 2D array per image)
% - can read a stack of images (returns 3D array)
%
% images are read according to x-y convention, use 'permute' flag to use
% Matlab y-x convention

% Thomas Deneux
% Copyright 2004-2012


% Input
if nargin<1
    fname = fn_getfile;
end
fname = cellstr(fname);
nimages = length(fname);

% first image
a = double(imread(fname{1})); 
if nargin<2
    a = permute(a,[2 1 3]); % Matlab (y,x) convention -> convention (x,y)
else
    if ~strcmp(flag,'permute'), error argument, end
end
if size(a,3)==3 && ~any(any(any(diff(a,1,3))))
    bw = true;
    a = a(:,:,1);
else
    bw = false;
end

% stack
if nimages>1
    if bw
        a(1,1,nimages) = 0;
    else
        a(1,1,1,nimages) = 0;
    end
    for i=2:nimages
        b = double(imread(fname{i}));
        if bw
            a(:,:,i) = b(:,:,1)';
        else
            a(:,:,:,i) = permute(b,[2 1 3]);
        end
    end
end

% make color image btw 0 and 1
if ~bw
    switch class(a)
        case 'uint8'
            a = double(a)/255;
        case 'double'
            nbyte = ceil(log2(max(a(:)))/8);
            switch nbyte
                case 1
                    a = a/255;
                case 2
                    a = a/65535;
                otherwise
                    disp('please help me')
                    keyboard
            end
        otherwise
            disp('please help me')
            keyboard
    end
end
    