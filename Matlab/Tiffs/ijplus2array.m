function array = ijplus2array(imp,newtype)
%IJPLUS2ARRAY Converts ImageJ ImagePlus object to MATLAB array
% array = ijplus2array(imp,newtype)
%
% by Vincent Bonin based on original code by Aaron Kerlin,

% changelog
% 01/25/08 vb corrects automatic unsigned integer cast
%          vb optimized speed and memory usage
%          vb allow bit depth conversion
% 03/25/10 vb optimized, 35% faster

% first convert to user-specified bit depth

% dummy object to disable automatic scaling in StackConverter
dummy = ij.process.ImageConverter(imp);
dummy.setDoScaling(0); % this is a static property

if nargin > 1
    if imp.getImageStackSize > 1
        converter = ij.process.StackConverter(imp); % doess not work for single images
    else
        converter = ij.process.ImageConverter(imp);
    end
    
    switch newtype
        case 'uint8'
            converter.convertToGray8;
        case 'uint16'
            converter.convertToGray16;  
        case 'single'
            converter.convertToGray32;        
        otherwise
            error('type not supported by imagej');
    end
end

% figure out 
switch imp.getBitDepth
    case 8
        arraytype = 'int8';
        srctype = 'uint8';
    case 16
        arraytype = 'int16';
        srctype = 'uint16';
    case 32
        arraytype = 'single'; % formerly 'single'
        srctype = 'single';
    otherwise
        error('unknown ij image type');
end

stack=imp.getImageStack();
nslices=stack.getSize();
width=stack.getWidth();
height=stack.getHeight();
pixelarray = stack.getImageArray();

%% VB optimized on 03/25/2010, 33 percent faster
this = zeros(width,height,srctype); %  permuted row and columns
array = zeros(height,width,nslices,srctype);
for islice = 1:nslices
    % because java has not signed types, matlab casts unsigned ij integers to signed
    this(:) = typecast(pixelarray(islice),srctype);
    array(:,:,islice)= this';
end;
array = reshape(array,[height,width,nslices]);

%% code prior to 03/25/2010
% array = zeros(width,height,nslices,srctype); % temporary type
% for islice = 1:nslices
%     % because java has not signed types, matlab casts unsigned ij integers to signed
%     array(:,:,islice)=reshape(typecast(pixelarray(islice),srctype),width,height);
% end
% 
% array=permute(array,[2 1 3]);
 
%%  original code by AK
%cellarray=cell(pixelarray); % converts java array into MATLAB cell array
%array=cell2mat(cellarray(1:nslices)); 

return;
