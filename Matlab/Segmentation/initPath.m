function [AllX, AllY, AllZ,AllDist] = initPath(Table2,Options)
%  [CenterX, CenterY, CenterZ] = initPath(Table,options)
% inputs: Table - containing X Y Z of the ROIs
%         Options - a struct with the following fields regarding the
%         referance stack and wanted ROIs ( = ? is default values):
%   1. Zstep = 1.5; %micrometer step in Z
%   2. UM_1X = 515; %micrometers at 1x zoom
%   3. Zoom = 2;
%   4. Pixels = 512;
%   5. InterpPoints = 3000; % number pf points to interpulate
%   6. LowpassWindow = 291 % window size to smooth the trace
%   7. xSize = 7 % size of FOV in um for X (assuming Y is the same)
%   8. zSize = 1.5 % size of FOV in um for X
%   9. StepSize = 7; % distance to place ROIs in um
%% parsing inputs

if nargin < 2
    throw(MException('initPath:Input','Two parameters needed'))
end
Table = Table2;
if isfield(Options,'Zstep')
    Zstep = Options.Zstep;
else
    Zstep = 1.5;
end
if isfield(Options,'UM_1X')
    UM_1X = Options.UM_1X;
else
    UM_1X = 515;
end
if isfield(Options,'Zoom')
    Zoom = Options.Zoom;
else
    Zoom = 2;
end
if isfield(Options,'RefPixels')
    RefPixels = Options.RefPixels;
else
    RefPixels = 1024;
end
if isfield(Options,'InterpPointsDist')
    InterpPointsDist = Options.InterpPointsDist;
else
    InterpPointsDist = 0.25;
end
if isfield(Options,'LowpassWindowUm')
    LowpassWindowUm = Options.LowpassWindowUm;
else
    LowpassWindowUm = 2;
end
if isfield(Options,'xSize')
    xSize = Options.xSize;
else
    xSize = 7;
end
if isfield(Options,'zSize')
    zSize = Options.zSize;
else
    zSize = 7;
end
if isfield(Options,'StepSize')
    StepSize = Options.StepSize;
else
    StepSize = 7;
end
if isfield(Options,'Plot')
    Plot = Options.Plot;
else
    Plot = 0;
end
if isfield(Options,'PadZ')
    PadZ = Options.PadZ;
else
    PadZ = 1;
end
%% Converting pixels to um
UM_Zoom = UM_1X./Zoom;
Pixel_Size = UM_Zoom./RefPixels;
Table.x = Table.x.*Pixel_Size;
Table.y = Table.y.*Pixel_Size;
Table.z = Table.z.*Zstep;
%% dealying with branches
AllDist = 0;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Find base branch
Base = find(Table.pNum == -1);
if isempty(Base) || length(Base) > 1
    throw(MException('initPath:base','Error with the tree: check disconnected branches'));
end
%% Go over the tree
AllX = [];
AllY = [];
AllZ = [];
while ~isempty(Base)
    %% get the base position
    X1 = Table.x(Base(1));
    Y1 = Table.y(Base(1));
    Z1 = Table.z(Base(1));
    %% get the parents
    Parents = find(Table.pNum == Base(1));
    Base(1) = [];
    if isempty(Parents)
        throw(MException('initPath:base','Error with the tree: Only one node!'));
    elseif length(Parents) > 1
        Base = Parents(2:end);
    end
    %% go over the rest of the branch
    Current = Parents(1);
    Flag = 1;
    while Flag
        X1 = [X1; Table.x(Current)];
        Y1 = [Y1; Table.y(Current)];
        Z1 = [Z1; Table.z(Current)];
        Parents = find(Table.pNum == Current);
        %         Table(Current,:) = [];
        if isempty(Parents)
            Flag = 0;
        elseif length(Parents) > 1
            Base = [Base; Parents(2:end)];
            Current = Parents(1);
        else
            Current = Parents(1);
        end
    end
    %% interpolation and smoothing
    cumDist=sum(sqrt(diff(X1).^2+diff(Y1).^2+ diff(Z1).^2));
    InterpPoints = round(cumDist./InterpPointsDist);
    LowpassWindow = round(cumDist./LowpassWindowUm);
    LowpassWindow = 2*floor(LowpassWindow/2)+1;
    V=interparc(InterpPoints,X1,Y1,Z1,'linear');
    X=sgolayfilt(V(:,1),1,LowpassWindow);
    Y=sgolayfilt(V(:,2),1,LowpassWindow);
    Z=sgolayfilt(V(:,3),1,LowpassWindow);
    %% distence function in 3d assuming square X=Y FOV
    Zprime = Z.*xSize./zSize; % weighing Z to ensure spacing with vertical dendrites
    dist=sqrt(diff(X).^2+diff(Y).^2+ diff(Zprime).^2);
    %% find center of volumes
    cDist = cumsum(dist);
    Vols = rem(cDist+StepSize/2,StepSize);
    [~,VolLocs] = findpeaks(Vols);
    CenterX = X(VolLocs);
    CenterY = Y(VolLocs);
    CenterZ = Z(VolLocs);
    %% Add scan fields by computing the local angle in Z
    if PadZ
        [CenterX,CenterY,CenterZ] = PadZplanes(CenterX,CenterY,CenterZ, Options);
    end
    %% add to global
    AllX = [AllX; CenterX];
    AllY = [AllY; CenterY];
    AllZ = [AllZ; CenterZ];
    AllDist = AllDist +cumDist;
end
%% Plot
if Plot
    figure
    scatter3(AllX,AllY,AllZ,'filled')
end