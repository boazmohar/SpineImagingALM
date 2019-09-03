function RegInfo=dendSeg3d(RegInfo,cell_vol_ups, ~)



if isfield(RegInfo.Cell,'labelimg_up')
    labelimg_up = RegInfo.Cell.labelimg_up;
    labelimg_up = imCellEditInteractive_3DMA_110407MA(cell_vol_up,(labelimg_up>0),[],1,1);
else
    labelimg_up = imCellEditInteractive_3DMA_110407MA(cell_vol_up,[],[],1,1);
end    
RegInfo.Cell.labelimg_up=labelimg_up;

labelimg_temp=labelimg_up((1:ceil(size(labelimg_up,1)./2))*2-1,(1:ceil(size(labelimg_up,2)./2))*2-1,(1:ceil(size(labelimg_up,3)./2))*2-1);
labelimg=zeros(size(labelimg_temp));
% reindex after removal of some sub-pixel cells
nCells_wrong=max(labelimg_temp(:));

counter=1;
for i=1:nCells_wrong
    cind=find(labelimg_temp==i);
    if ~isempty(cind)
        [x,y,z]=ind2sub(size(labelimg_temp),cind);
        for j=1:length(x)
            labelimg(x(j),y(j),z(j))=counter;
        end
        counter=counter+1;
    end
end

nCells=counter-1;
binarymask_up=(labelimg_up>0);
binarymask=binarymask_up([1:ceil(size(labelimg_up,1)./2)]*2-1,[1:ceil(size(labelimg_up,2)./2)]*2-1,[1:ceil(size(binarymask_up,3)./2)]*2-1);  
strr=regionprops(labelimg,'Centroid', 'Area');
centroidtmp={strr.Centroid};
Centroid = cell2mat(centroidtmp')';
Area = [strr.Area]; 

RegInfo.Cell.binarymask = binarymask;
RegInfo.Cell.labelimg = labelimg;
RegInfo.Cell.nCells = nCells;
RegInfo.Cell.centroid = Centroid;
RegInfo.Cell.area = Area;