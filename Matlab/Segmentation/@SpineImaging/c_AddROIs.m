function c_AddROIs(obj)
hSI = evalin('base','hSI');
hSI.hRoiManager.roiGroupMroi.clear();
hSI.hRoiManager.mroiEnable = true;
% hSI.hDisplay.displayVolumeLock = true;          %set display to show all z planes
hSI.hFastZ.enable = true;
hSI.hFastZ.numVolumes = obj.OptionsStruct.numVolumes;
% hSI.hFastZ.settlingTime = obj.OptionsStruct.settlingTime;
hSI.hFastZ.actuatorLag = obj.OptionsStruct.settlingTime;
hSI.hFastZ.volumePeriodAdjustment = -obj.OptionsStruct.volPeriodAdj/1000000;
% hSI.hFastZ.discardFlybackFrames = false;
hSI.hFastZ.flybackTime = 0;
hSI.hFastZ.waveformType = 'step';
hSI.hFastZ.useArbitraryZs = true;
userZs = obj.All.z';
while length(unique(userZs)) < length(userZs)
    u=unique(userZs);
    n=histc(userZs,u);
    problems = u(n>1);
    if ~isempty(problems)
        problems = problems(1);
        index = find(userZs==problems,2);
        if length(index) > 1
            userZs(index(2)) = userZs(index(2))+0.01;
        end
    end
end
hSI.hFastZ.userZs = userZs;  %Z to scan in
% setting frame flyto and flayback in ms

% hScan = hSI.hResScan
hScan = hSI.hScan2D;

hScan.flybackTimePerFrame = obj.OptionsStruct.settlingTime;
hScan.flytoTimePerScanfield = obj.OptionsStruct.settlingTime;
% setting the fill fraction for mROI
hScan.fillFractionSpatial = obj.OptionsStruct.mROIFillFactor;
%setting frames per file to be the number of ROIs
if isfield(obj.OptionsStruct,'DenseStackMode') && ...
          obj.OptionsStruct.DenseStackMode
    hScan.logFramesPerFile = length(obj.All.z);
    hSI.hFastZ.numVolumes = 200;
else
    hScan.logFramesPerFile = Inf;
end
% change beam lead time and power
%
% hScan.beamClockDelay = obj.OptionsStruct.BeamLeadTime;
hScan.beamClockDelay = obj.OptionsStruct.BeamLeadTime/10^6;
% hScan.beamClockExtend set to 1us or more

hSI.hBeams.powers = obj.OptionsStruct.Power;
if  ~isempty(obj.OptionsStruct.ZPowerAdjLength)
    hSI.hBeams.pzAdjust = 1;
    hSI.hBeams.lengthConstants = obj.OptionsStruct.ZPowerAdjLength;
else
     hSI.hBeams.pzAdjust = 0;
end
%% check for low power fields
Branches = floor(obj.OptionsStruct.LowBranches);
Id = [];
for i = 1:length(Branches)
    txt = num2str(obj.OptionsStruct.LowBranches(i));
    dot = strfind(txt, '.');
    if ~isempty(dot) && isempty(obj.OptionsStruct.MaxBranchDist)
        throw(MException('SpineImaging:c_SelectBranch', ...
            'Do not use dot notation if MaxBranchDist is off'))
    end
    if ~isempty(dot)
        txt = txt(dot+1:end);
        SubBranch = zeros(length(txt),1);
        for s = 1:length(txt)
            Id = [Id Branches(i) + str2double(txt(s))/10];
        end
    else
        Id = [Id Branches(i)];
    end
end
%% Create the roi group
myGroup = scanimage.mroi.RoiGroup(obj.OptionsStruct.GroupName);
% dendrites
dendriteRoi = scanimage.mroi.Roi;
dendriteRoi.discretePlaneMode = true;
%low power dendrites (cell body)
lowRoi = scanimage.mroi.Roi;
lowRoi.discretePlaneMode = true;
lowRoi.powers = obj.OptionsStruct.Power;
% blanked flyto and flyback fields
blankedRoi = scanimage.mroi.Roi;
blankedRoi.discretePlaneMode = true;
blankedRoi.display = false;
blankedRoi.powers = 0;
for i = 1:length(obj.Roi.x)
    % sf1 should be: 1) 4 elemnt vector for size in relative units
                    %2) rotation angle (only 0 for resonant)
                    %3) number of pixels and lines
    sf1 = scanimage.mroi.scanfield.fields.RotatedRectangle(...
        [obj.Roi.x(i) obj.Roi.y(i) obj.Roi.xSize obj.Roi.ySize],0,...    
        [obj.OptionsStruct.ImagingPixels,obj.OptionsStruct.ImagingLines]);
    if obj.TrueFieldsMask(i)                    % a real field
        if ~isempty(Id) && sum(obj.All.Id(i) == Id)
            lowRoi.add(userZs(i),sf1);                  % low power field
        else
            dendriteRoi.add(userZs(i),sf1);             % regular field
        end
    else
         blankedRoi.add(userZs(i),sf1);         % a blanked field
    end
end
myGroup.add(dendriteRoi);
myGroup.add(blankedRoi);
myGroup.add(lowRoi);
hSI.hRoiManager.roiGroupMroi = myGroup;
hSI.hDisplay.resetActiveDisplayFigs();
evalin('base','hSICtl.editImagingRoiGroup');
if ~obj.OptionsStruct.ShowScanField
    hSI.hDisplay.roiDisplayEdgeColor='none';
    hSI.hDisplay.roiProjectionDisplayEdgeColor='none';
else
    hSI.hDisplay.roiDisplayEdgeColor='b';
    hSI.hDisplay.roiProjectionDisplayEdgeColor='b';
end
