function [CenterX, CenterY,CenterZ] = c_PadZ(obj,CenterX,CenterY,CenterZ)
% [CenterX,CenterY,CenterZ] = SegmentDendrite(CenterX,CenterY,CenterZ, Options)
% Add scan fields by computing the local angle in Z
% % Inputs: Centers of X Y Z of the ROIs
%         Options - a struct with the following fields regarding the
%         referance stack and wanted ROIs ( = ? is default values):
%   1. Thres1 threshold above wich a plane is added in degress
%   2. Thres2 same - adds 2 planes +-1.5
%   3. Plot - to display figures
% Outputs: Tables each containg centers of seperate dendrites
%% compute angle
dx = diff(CenterX);
dy = diff(CenterY);
dz = diff(CenterZ);
Angle = abs(atand(dz./(sqrt(dx.^2 + dy.^2))));
%% which ROIs to duplicate
Half = obj.OptionsStruct.zSize./2;
Full = obj.OptionsStruct.zSize;
AngleThres1 = obj.OptionsStruct.AngleThres1;
AngleThres2 = obj.OptionsStruct.AngleThres2;
Add1 = Angle>AngleThres1 & Angle <=AngleThres2;
Add2 = Angle<=AngleThres1;
if ~isempty(Add1) && ~isempty(Add2)
    Add1 = [Add1; Add1(end)];
    Add2 = [Add2; Add2(end)];
    CenterX2 = [CenterX;    CenterX(Add1);       CenterX(Add2);     CenterX(Add2)];
    CenterY2 = [CenterY;    CenterY(Add1);       CenterY(Add2);     CenterY(Add2)];
    CenterZ2 = [CenterZ;    CenterZ(Add1)+Half;  CenterZ(Add2)+Full; CenterZ(Add2)-Full];
    CenterZ2(Add1) = CenterZ2(Add1) - Half;
    
    %% return results
    CenterX = CenterX2;
    CenterY = CenterY2;
    CenterZ = CenterZ2;
end