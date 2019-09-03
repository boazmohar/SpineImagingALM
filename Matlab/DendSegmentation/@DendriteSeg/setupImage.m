function setupImage(obj)
%% get min max vals
volume = obj.display.volume;
IMG_CONTRAST_INV = 1 - obj.display.Im_contrast;
MINMAX0 = [min(volume(:)) max(volume(:))];
MINMAX = [(MINMAX0(1) + IMG_CONTRAST_INV*diff(MINMAX0)) ((MINMAX0(2) - IMG_CONTRAST_INV*diff(MINMAX0)))]; 
%% get size
[nRows, nCols, nZ] = size(volume);
obj.display.nRows = nRows;
obj.display.nCols = nCols;
obj.display.nZ = nZ;
obj.display.MINMAX = MINMAX;
%% dispaly first image
obj.display.brightset = 1;
% create an initial empty mask if one isn't loaded
obj.cells.maskSpine = false(nRows, nCols, nZ);
obj.cells.maskBorder = false(nRows, nCols, nZ);
obj.cells.labelBorder = zeros(nRows, nCols, nZ);
%%
plane = obj.display.volume(:,:,obj.display.Z);
scaledImage = obj.normImage(plane);
mask = obj.cells.maskSpine(:,:,obj.display.Z); 
sahdedImage = imShade(scaledImage,mask);
obj.handles.imH = imagesc(sahdedImage,MINMAX);
set(obj.handles.axH, 'DataAspectRatio', [1 1 1])
%% add pixel Info to track mouse location during clicking
pxInfoH = impixelinfoval(obj.handles.figH,obj.handles.imH);
set(pxInfoH,'position',[900,5,300,20])
obj.handles.pxInfoH = pxInfoH;
%% se
obj.display.se = strel('square',3);