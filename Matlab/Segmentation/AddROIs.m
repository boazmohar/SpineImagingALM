function AddROIs(UpperX, LeftY,CenterZ, RoiSize,Options)

if nargin < 4
    throw(MException('MATLAB:ambiguousSyntax','Four parameters needed'))
end
if isfield(Options,'numVolumes')
    numVolumes = Options.numVolumes;
else
   numVolumes =200;
end
if isfield(Options,'settlingTime')
    settlingTime = Options.settlingTime;
else
    settlingTime = 0.01;
end
if isfield(Options,'GroupName') 
    GroupName = Options.GroupName;
else
    GroupName = 'defult group name';
end
if isfield(Options,'Pixels') 
    Pixels = Options.linesPerField;
else
    Pixels = 30;
end
hSI = evalin('base','hSI');
hSI.hRoiManager.mroiEnable = true;
hSI.hDisplay.displayVolumeLock = true;          %set display to show all z planes
hSI.hFastZ.enable = true;
hSI.hFastZ.numVolumes = numVolumes;
hSI.hFastZ.settlingTime = settlingTime;
hSI.hFastZ.discardFlybackFrames = false; 
hSI.hFastZ.waveformType = 'step';
hSI.hFastZ.useArbitraryZs = true;              
hSI.hFastZ.userZs = CenterZ;  %Z to scan in 
%% Create the roi group
myGroup = scanimage.mroi.RoiGroup(GroupName);
dendriteRoi = scanimage.mroi.Roi;
for i = 1:length(UpperX)
    sf1 = scanimage.mroi.scanfield.fields.RotatedRectangle(...
    [UpperX(i) LeftY(i) RoiSize RoiSize],...    %Rect [x y width height]. Units are scan fraction (0-1) of the reference field of view
    0,...                                       %Rotation angle. Only zero supported for resonant scanning
    [Pixels,Pixels]);                           %Pixel resolution
    dendriteRoi.add(CenterZ(i),sf1);            %add scanfields add appropriate z planes.
    clear sf1;
end
myGroup.add(dendriteRoi);

%%
hSI.hRoiManager.roiGroupMroi = myGroup;
return
hSI.editImagingRoiGroup();   %open the roi group editor to view our handywork
%% clean up
clear ans myGroup dendriteRoi;