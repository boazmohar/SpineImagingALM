function [CenterX,CenterY,CenterZ] = PadZplanes(CenterX,CenterY,CenterZ, Options)
% [CenterX,CenterY,CenterZ] = SegmentDendrite(CenterX,CenterY,CenterZ, Options)
% Add scan fields by computing the local angle in Z
% % Inputs: Centers of X Y Z of the ROIs
%         Options - a struct with the following fields regarding the
%         referance stack and wanted ROIs ( = ? is default values):
%   1. Thres1 threshold above wich a plane is added in degress
%   2. Thres2 same - adds 2 planes +-1.5
%   3. Plot - to display figures
% Outputs: Tables each containg centers of seperate dendrites
%%
if nargin < 4
    throw(MException('MATLAB:ambiguousSyntax','4 parameters needed'))
end
if isfield(Options,'AngleThres1')
    AngleThres1 = Options.AngleThres1;
else
    AngleThres1 =22.5;
end
if isfield(Options,'AngleThres2')
    AngleThres2 = Options.AngleThres2;
else
    AngleThres2 = 67.5;
end
if isfield(Options,'Plot')
    Plot = Options.Plot;
else
    Plot = 0;
end
%% compute angle
dx = diff(CenterX);
dy = diff(CenterY);
dz = diff(CenterZ);
Angle = abs(atand(dz./(sqrt(dx.^2 + dy.^2))));
%% which ROIs to duplicate
Add1 = Angle>AngleThres1 & Angle <=AngleThres2;
Add2 = Angle<=AngleThres1;
Add1 = [Add1; Add1(end)];
Add2 = [Add2; Add2(end)];
CenterX2 = [CenterX;    CenterX(Add1);       CenterX(Add2);     CenterX(Add2)];
CenterY2 = [CenterY;    CenterY(Add1);       CenterY(Add2);     CenterY(Add2)];
CenterZ2 = [CenterZ;    CenterZ(Add1)+0.75;  CenterZ(Add2)+1.5; CenterZ(Add2)-1.5];
CenterZ2(Add1) = CenterZ2(Add1) -0.75;
%% ploting
if Plot
    figure;
    clf
    scatter3(CenterX,CenterY,CenterZ,'filled')
    hold on
    scatter3(CenterX2,CenterY2,CenterZ2,'filled')
end
%% return results
CenterX = CenterX2;
CenterY = CenterY2;
CenterZ = CenterZ2;