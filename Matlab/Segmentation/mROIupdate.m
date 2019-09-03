function mROIupdate(Options)

global CenterX CenterY CenterZ
if Options.ReadFromFile
    Table = readSWC(Options);
    [CenterX, CenterY, CenterZ,AllDist] = initPath(Table,Options);
else
    AllDist = -1;
    CenterX(isnan(CenterZ)) = [];
    CenterY(isnan(CenterZ)) = [];
    CenterZ(isnan(CenterZ)) = [];
    if Options.PadZ
        [CenterX,CenterY,CenterZ] = PadZplanes(CenterX,CenterY,CenterZ, Options)
    end
end
    
%% Get an initial path
axes(evalin('base','aHandle'));
scatter3(CenterX, CenterY, CenterZ,'filled')
linkdata on;
axis equal;
%%
framerate = 1/(Options.linesPerField/16000*length(CenterX));
title(['Total length: ' num2str(AllDist,3) 'um, Vol rate: ' num2str(framerate,3) 'Hz']);