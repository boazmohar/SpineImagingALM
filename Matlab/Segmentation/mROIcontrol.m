close all; clear all; clc

global CenterZ CenterY CenterX;
fHandle = figure(1);
fHandle.Units = 'Normalized';
fHandle.Position = [0.05 0.05 0.9 0.85];
pHandle = uipanel('parent',fHandle','position',[0 0 0.3 1]);
aHandle = axes('parent',fHandle,'position',[0.32 0.05 0.625 0.9]);
fHandle.Units = 'pixels';
%
OptionsStruc.filename = 'Z:\Analysis_share\segmentation sample\WRBM6\Cell1_3dent_v1.swc';
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
OptionsStruc.Plot = 0;               % whether to plot
Spec.Plot = 'logical';
OptionsStruc.PadZ = 1;            % If to pad Z with extra planes
Spec.PadZ = 'logical';
OptionsStruc.AngleThres1 = 45/2;          % threshold to add 1 more plane (seperated by 1.5um)
OptionsStruc.AngleThres2 = 45/2 + 45;     % threshold to add 2 more plane (seperated by 1.5um)
OptionsStruc.numVolumes = 200;       % How many volumes to do ?
OptionsStruc.settlingTime = 0.000126;    % seteling time
OptionsStruc.GroupName = 'Cell x from WRxx'; %mROI group name
OptionsStruc.Pixels = 22;            % number of pixels and lines of the ROI
Spec.Pixels = 'slider 1 128 1 %d';
OptionsStruc.maxAccel = 0.5;
OptionsStruc.maxVel = 18;
OptionsStruc.lineRate = 16000;
OptionsStruc.linesPerField = 24;
OptionsStruc.targetAccuracy = 0.75;
OptionsStruc.zOptIter = 10;
OptionsStruc.volPeriodAdj = 500; %in microseconds
OptionsStruc.ReadFromFile = 1; %in microseconds
Spec.ReadFromFile = 'logical';

OptionsObj = fn_control(OptionsStruc,@mROIupdate,Spec,pHandle); % to change a value change the GUI
% Pos = OptionsObj.hp.Position;
% set(OptionsObj.hp,'Position', [10 10 Pos(3) Pos(4)]);
mROIupdate(OptionsObj.s);