function [UpperX, LeftY, RoiSize] = ConvertROIs(CenterX, CenterY, Options)
%  UpperX, LeftY, RoiSize] = ConvertROIs(CenterX, CenterY, Options)
% inputs: The center of the ROIs in X and Y
%         Options - a struct with the following fields regarding the
%         referance stack and wanted ROIs ( = ? is default values):
%   1. UM_1X = 515; %micrometers at 1x zoom
%   2. Zoom = 2;
%   3. xSize = 7 % size of FOV in um for X (assuming Y is the same)
%   4. zSize = 1.5 % size of FOV in um for X
% Outputs: The upper (X) left (Y) corners of the ROI and it's size in
% relative units to the 1x scan field.
if nargin < 3
    throw(MException('MATLAB:ambiguousSyntax','Three parameters needed'))
end
if isfield(Options,'Zoom')
    Zoom = Options.Zoom;
else
   Zoom =2;
end
if isfield(Options,'UM_1X')
    UM_1X = Options.UM_1X;
else
    UM_1X = 515;
end
if isfield(Options,'xSize')
    xSize = Options.xSize;
else
    xSize = 7;
end
if isfield(Options,'ySize')
    ySize = Options.ySize;
else
    ySize = 7;
end
ZoomOffset = (1 - 1/Zoom)/2;
CenterXNorm = CenterX./UM_1X + ZoomOffset;  % Normalizing & shiffting becuse of zoom
CenterYNorm = CenterY./UM_1X + ZoomOffset; 
UpperX = CenterXNorm - xSize./2./UM_1X;     % Adding the ROI size
LeftY = CenterYNorm - ySize./2./UM_1X;
RoiSize =  xSize./UM_1X;
if sum(UpperX < 0)
    disp('Warning: negative ROI x value found')
    UpperX(UpperX<0) = 0;
end

if sum(LeftY < 0)
    disp('Warning: negative ROI y value found')
    LeftY(LeftY<0) = 0;
end