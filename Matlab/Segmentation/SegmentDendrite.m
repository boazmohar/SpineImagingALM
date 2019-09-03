function Tables = SegmentDendrite(Table, Options)
%  Tables = SegmentDendrite(Table, Options) segments the dendrite if it has
%  any branch points into seperate trajectories for InitPath
% Inputs: Table - containing X Y Z of the ROIs
%         Options - a struct with the following fields regarding the
%         referance stack and wanted ROIs ( = ? is default values):
%   1. Zstep = 1.5; %micrometer step in Z
% Outputs: Tables each containg centers of seperate dendrites
%% input handeling
if nargin < 2
    throw(MException('MATLAB:ambiguousSyntax','Two parameters needed'))
end
if isfield(Options,'Zoom')
    Zoom = Options.Zoom;
else
   Zoom =2;
end
%%
