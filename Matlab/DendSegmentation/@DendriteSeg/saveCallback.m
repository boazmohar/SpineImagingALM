function saveCallback(obj,figH,eventData)
set(obj.handles.ustatH,'String','Saving')
pause(0.1)
%% reindex after removal of some sub-pixel cells
labelimg_up = obj.cells.labelSpine;
labelimg_up(obj.cells.maskBorder) = 0;
labelimg_temp=labelimg_up((1:ceil(size(labelimg_up,1)./2))*2-1,(1:ceil(size(labelimg_up,2)./2))*2-1,(1:ceil(size(labelimg_up,3)./2))*2-1);
labelimg=zeros(size(labelimg_temp));
nCells_wrong=max(labelimg_temp(:));
goodCells = [];
counter=1;
for i=1:nCells_wrong
    cind=find(labelimg_temp==i);
    if ~isempty(cind)
        [x,y,z]=ind2sub(size(labelimg_temp),cind);
        for j=1:length(x)
            labelimg(x(j),y(j),z(j))=counter;
        end
        goodCells = [goodCells i];
        counter=counter+1;
    end
end
%% calculating
nCells=counter-1;
binarymask_up=(labelimg_up>0);
binarymask=binarymask_up([1:ceil(size(labelimg_up,1)./2)]*2-1,[1:ceil(size(labelimg_up,2)./2)]*2-1,[1:ceil(size(binarymask_up,3)./2)]*2-1);
strr=regionprops(labelimg,'Centroid', 'Area');
centroidtmp={strr.Centroid};
Centroid = cell2mat(centroidtmp')';
Area = [strr.Area];
%% saveing to struct
cells = obj.cells;
cells.downsampled.nCells = nCells;
cells.downsampled.binarymask = binarymask;
cells.downsampled.Centroid = Centroid;
cells.downsampled.Area = Area;
cells.downsampled.labelimg = labelimg;
%% saving to disk
[FileName,PathName] = uiputfile();
% ('cells','mask.mat')
%%
newCell = struct();
newCell.labelimg = uint16(cells.downsampled.labelimg);
newCell.dend = cells.branches;
all = [cells.all{goodCells}];
newCell.dendNum = [all(:).dendNum];
newCell.quality = {all(:).Quality};
newCell.note = {all(:).note};
newCell.dendTable = table2array(cells.dendriteTable);

%%
save([PathName FileName '_python.mat'],'-struct','newCell','-v6')
save([PathName FileName '.mat'],'cells','-v7.3')
%%
obj.saving.changed = 0;
set(obj.handles.ustatH,'String','Ready')
end