function c_SetDefaultOptions(obj)
%% ref image properties
OptionsStruc.filename = '';                 % SWC filename to load
Spec.filename = 'file';
OptionsStruc.Stack = '';                    % Stack filename to load
Spec.Stack = 'file';
OptionsStruc.typeNumbers = -1;              % What segment type to get from SWC file
OptionsStruc.Zstep = 1.6;                   % micrometer step in Z
OptionsStruc.UM_1X = 535;                   % micrometers at 1x zoom
OptionsStruc.Zoom = 1;                      % zoom for referance stack
OptionsStruc.RefPixels = 1024;              % Number of pixels for referace stack
OptionsStruc.RefFillFactor = 0.891;         % Spatial fill factor of reference image
%% scanfield properties
OptionsStruc.InterpPointsDist = 0.25;       % um distance between points to inturpalate
OptionsStruc.LowpassWindowUm = 3;           % window size to smooth the trace
OptionsStruc.xSize = 24;                    % size of FOV in um for X (assuming Y is the same)
OptionsStruc.ySize = 12;      
OptionsStruc.zSize = 3;                     % size of FOV in um for Z
OptionsStruc.zAspectRatio = 3;            % used in determning distence function 
OptionsStruc.StepSize = 12;                 % distance to place ROIs in um
OptionsStruc.AngleThres1 = 45/2;            % threshold to add 1 more plane (seperated by 1.5um)
OptionsStruc.AngleThres2 = 45/2 + 45;       % threshold to add 2 more plane (seperated by 1.5um)
OptionsStruc.ImagingPixels = 72;           % number of pixels and lines of the ROI
OptionsStruc.ImagingLines = 36;             % number of pixels and lines of the ROI
OptionsStruc.FlyLines = 4;                  % fly back lines
OptionsStruc.BeamLeadTime = 0.8;            % beam lead time in mROI is different
% Spec.Pixels = 'slider 1 128 1 %d';
%% Z optimization
OptionsStruc.maxAccel = 20;                  % Max acceleration of Z
OptionsStruc.maxVel = 40;                   % Max velocity of Z
OptionsStruc.lineRate = 16000;              % x resonant line rate
OptionsStruc.offset = 0;                    % offset in microns
OptionsStruc.targetAccuracy = 0.75;         % max diff in z (um) from trajectory
OptionsStruc.zOptRep = 5;                  % Repetitions to average
OptionsStruc.zOptIter = 10;                 % Max iterations
OptionsStruc.volPeriodAdj = 2000;            % in microseconds
OptionsStruc.settlingTime = 0.000252;       % seteling time
OptionsStruc.ZmaxTravel = 1600;              % max travel of the z axis
OptionsStruc.mROIFillFactor = 0.891;        % mROI fill factor
%% Ploting
OptionsStruc.PlotReconstraction =1;         % plot table centers
Spec.PlotReconstraction = 'logical';
OptionsStruc.PlotCenters =1;                % plot scanfield centers
Spec.PlotCenters = 'logical';
OptionsStruc.PlotLine =1;                   % plot interpolated center line
Spec.PlotLine = 'logical';
OptionsStruc.ShowScanField = 1;             % Show Scan Field broders
Spec.ShowScanField = 'logical';
OptionsStruc.ShowZoptoFields = 1;             % Show Scan Field added in zOpto
Spec.ShowZoptoFields = 'logical';
%% processing stages
OptionsStruc.ZoptDampen= true;             % to add control points to Z optimization
OptionsStruc.PadZ = false;                  % If to pad Z with extra planes
OptionsStruc.SelectedBranches = [];         % which branches to scan
Spec.SelectedBranches = 'xdouble 5 [1]'; 
OptionsStruc.LowBranches = [];              % which branches to scan low sampleing
Spec.LowBranches = 'xdouble 5 [1]'; 
OptionsStruc.MaxBranchDist = [];            % max dist a branch could be befor sub deviding
Spec.MaxBranchDist = 'xdouble 1 [50]'; 
% OptionsStruc.BrushingMode = false;          % Reset to centers from file
OptionsStruc.ZOptimize = false;             % Preform Z optimization
OptionsStruc.Power = 15;                    % Power at 0um
OptionsStruc.ZPowerAdjLength = [];          % Inct
Spec.ZPowerAdjLength = 'xdouble 5 [200]'; 
OptionsStruc.DenseStackMode = false;        % chages mode to dense stack
%% Other
OptionsStruc.numVolumes = 100000 ;              % How many volumes to do 
OptionsStruc.GroupName = 'Cell x';% mROI group name

%% Add to the object
obj.OptionsStruct = OptionsStruc;
obj.OptionsSpec = Spec;