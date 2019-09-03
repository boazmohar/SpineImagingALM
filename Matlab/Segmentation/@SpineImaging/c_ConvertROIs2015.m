function c_ConvertROIs2015(obj)
% 
% Outputs: The upper (X) left (Y) corners of the ROI and it's size in
% relative units to the 1x scan field.
%% normalizeing by zoom
ZoomOffset = (1 - 1/obj.OptionsStruct.Zoom)/2;
CenterXNorm = obj.All.x./obj.OptionsStruct.UM_1X + ZoomOffset;  % Normalizing & shiffting becuse of zoom
CenterYNorm = obj.All.y./obj.OptionsStruct.UM_1X + ZoomOffset; 
%% Normalizng by fill factor
FillFactorOffset = (1 - obj.OptionsStruct.RefFillFactor)/2;
CenterXNorm = CenterXNorm*obj.OptionsStruct.RefFillFactor + FillFactorOffset;
CenterYNorm = CenterYNorm*obj.OptionsStruct.RefFillFactor + FillFactorOffset;
%% Converting to left upper corner
obj.Roi.x = CenterXNorm - obj.OptionsStruct.xSize./2./obj.OptionsStruct.UM_1X; 
obj.Roi.y = CenterYNorm - obj.OptionsStruct.ySize./2./obj.OptionsStruct.UM_1X;
obj.Roi.ySize =  obj.OptionsStruct.ySize./obj.OptionsStruct.UM_1X;
obj.Roi.xSize =  obj.OptionsStruct.xSize./obj.OptionsStruct.UM_1X;
if sum(obj.Roi.x < 0)
    disp('Warning: c_ConvertROIs - negative values in x')
    obj.Roi.x(obj.Roi.x<0) = 0;
end
if sum(obj.Roi.y < 0)
    disp('Warning: c_ConvertROIs - negative values in y')
    obj.Roi.y(obj.Roi.y<0) = 0;
end