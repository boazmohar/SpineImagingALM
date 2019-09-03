function imp = ijarray2plus(array, typestr, ijName)
%IJARRAY2PLUS Converts MATLAB array to ImageJ ImagePlus object
%   imp = ijarray2plus(array, typestr)
%   array size [nRows,nCols,nFrames,nPlanes]
% $Id: ijarray2plus.m 475 2009-03-20 20:43:13Z kenichi $

if nargin < 3, 
    % get caller's name for array
    ijName = inputname(1);
    if isempty(ijName), ijName = 'stack'; end
end
    
%%%

[height,width,nframes,nplanes]=size(array);

arraytype = class(array);
% don't convert by default
if nargin < 2 || isempty(typestr), typestr = arraytype; end  
    

%% Handle different bit depths and RGB images
if nplanes > 1
    assert(nplanes == 3, 'color channel (4th dim) must be size 3 -- RGB');

    assert(strcmp(arraytype, 'uint8'), 'RGB arrays must be 8 bit for ImageJ');
    
    process = @ij.process.ColorProcessor;
else
    switch arraytype
      case {'uint8'}
        process = @ij.process.ByteProcessor;
      case {'uint16'}
        process = @ij.process.ShortProcessor;
      case {'uint32','single','double'}
        process = @ij.process.FloatProcessor;        
      otherwise
        error('type note supported');
    end

    if any(strcmp(arraytype,{'uint32','double'}))
        array = single(array); % uint32, double not supported in imagej
        arraytype = 'single';
    end
end

%% make a stack
stack = ij.ImageStack(width,height);

%% go through each frame and convert
for i=1:nframes
    ip = process(width,height);
    if nplanes == 3
        %% rgb
        rPix = reshape(array(:,:,i,1)',width*height,1);
        gPix = reshape(array(:,:,i,2)',width*height,1);        
        bPix = reshape(array(:,:,i,3)',width*height,1);                
        ip.setRGB(rPix, gPix, bPix);
    else
        %% one color
        pixels = reshape(array(:,:,i)',width*height,1);
        ip.setPixels(pixels);                   
    end
    stack.addSlice('frame', ip);
end

clear array;

imp = ij.ImagePlus(ijName,stack);

if nplanes == 3
    %% note: RGB stack does not need conversion.
    %% Do the other bit depths need it too?  I don't know.  Feel free to
    %% remove if its safe to do it - MH 080730
    return
end


% dummy object to disable automatic scaling in StackConverter
dummy = ij.process.ImageConverter(imp);
dummy.setDoScaling(0); % this is a static property

if ~strcmp(arraytype,typestr)
    if nframes > 1
        converter = ij.process.StackConverter(imp);
    else
        converter = ij.process.ImageConverter(imp);
    end
    switch typestr
        case 'uint8';
            converter.convertToGray8;
        case 'uint16';
            converter.convertToGray16;
        case {'single','uint32','double'};
            converter.convertToGray32;
    end
end

return;
