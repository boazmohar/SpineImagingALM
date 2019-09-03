function swcCallback(obj,figH,eventData)
%% read swc file
obj.dendReadSWC();
%% detect branches
obj.dendAddBranches();
%% recalculate mask
obj.dendMaskUpdate();