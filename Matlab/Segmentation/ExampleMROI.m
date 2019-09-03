%% Work flow example for mROI
%clear all; clc; close all;
%% setting options struct
OptionsStruc.filename = 'T:\Analysis_share\segmentation sample\WRBM6\Cell1_3dent_v1.swc';
Spec.filename = 'file';
OptionsStruc.typeNumbers = -1;
OptionsStruc.Zstep = 1.5;            % micrometer step in Z
OptionsStruc.UM_1X = 515;            % micrometers at 1x zoom
OptionsStruc.Zoom = 1;               % zoom for referance stack
OptionsStruc.RefPixels = 1024;       % Number of pixels for referace stack
OptionsStruc.InterpPointsDist = 0.25;% um distance between points to inturpalate
OptionsStruc.LowpassWindowUm = 2;    % window size to smooth the trace
OptionsStruc.xSize = 7;              % size of FOV in um for X (assuming Y is the same)
OptionsStruc.zSize = 1.5;            % size of FOV in um for X
OptionsStruc.StepSize = 7;           % distance to place ROIs in um
OptionsStruc.Plot = 1;               % whether to plot
Spec.Plot = 'logical';
OptionsStruc.PadZ = 1;            % If to pad Z with extra planes
Spec.PadZ = 'logical';
OptionsStruc.AngleThres1 = 45/2;          % threshold to add 1 more plane (seperated by 1.5um)
OptionsStruc.AngleThres2 = 45/2 + 45;     % threshold to add 2 more plane (seperated by 1.5um)
OptionsStruc.numVolumes = 200;       % How many volumes to do ?
OptionsStruc.settlingTime = 0.000126;    % seteling time
OptionsStruc.GroupName = 'Cell x from WRxx'; %mROI group name
OptionsStruc.Pixels = 22;            % number of pixels and lines of the ROI
OptionsStruc.maxAccel = 0.5;
OptionsStruc.maxVel = 18;
OptionsStruc.lineRate = 16000;
OptionsStruc.linesPerField = 24;
OptionsStruc.targetAccuracy = 0.75;
OptionsStruc.zOptIter = 10;
OptionsStruc.volPeriodAdj = 500; %in microseconds
Spec.Pixels = 'slider 1 128 1 %d';
OptionsObj = fn_control(OptionsStruc,Spec); % to change a value change the GUI
%% read data from SWC file
Table = readSWC(OptionsObj.s);
%% Get an initial path
[CenterX, CenterY, CenterZ,AllDist] = initPath(Table, OptionsObj.s);
%% Z Re-profile
clear FieldCenters
FieldCenters(:,1)=CenterX;
FieldCenters(:,2)=CenterY;
FieldCenters(:,3)=CenterZ;
[Zwave,NewFieldCenters,TrueFieldsMask,TrueFieldsCenterSamp]=ZOptimize(FieldCenters,OptionsObj.s);
NewCenterX=NewFieldCenters(:,1);
NewCenterY=NewFieldCenters(:,2);
NewCenterZ=NewFieldCenters(:,3);
% input number of lines, Center of X Y Z scan fields, max accelertion vector
% returns Z profile and new set of Centers for X Y Z with padding if needed
% make sure we don't have two exect same Z positions
%% Convert Centers to SI mROI = [x y width height]. Units are scan fraction (0-1) of the reference field of view
[UpperX, LeftY, RoiSize] = ConvertROIs(NewCenterX, NewCenterY, OptionsObj.s); % check for edge!!!
%% Optimizting Z piezo response to hit target positions 
[CmdWave,maxErrorMicron]=PiezoOptIterative(Zwave,TrueFieldsCenterSamp,OptionsObj.s);
maxError
%% Getting SI ready and adding ROIs (add Z waveform setup)sc
    AddROIs(UpperX, LeftY, NewCenterZ, RoiSize, OptionsObj.s);
