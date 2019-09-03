function addCell3d(obj)
%% get prop from obj
obj.handles.ustatH.String ='Adding mask';
Z1 = obj.display.Z;
nZ = obj.display.nZ;
fovImg_vol = obj.display.volume;
bwCurr = squeeze(obj.cells.maskSpine(:,:,Z1));
fovImg = squeeze(obj.display.volume(:,:,Z1));
bwNew_vol0 = false(size(obj.cells.maskSpine));
bwNew_vol = obj.cells.maskSpine;
diskR = obj.current.diskR;
%% fill in whole volume, use same mean as mean from current plane to get whole cell:
[~,~,cellMean] = obj.subAddCell(bwCurr, fovImg);
%% look for other Z planes 
for step_Z = [1 -1] % 1 --> positive Z direction, -1 negative Z direction
    STOP = 0;
    if step_Z == 1
        Z1_use = Z1;
    elseif step_Z == -1
        Z1_use = Z1 - 1;
    end
    while STOP ~= 1 && (Z1_use > 0) && (Z1_use <= nZ)
        fovImg_USE = squeeze(fovImg_vol(:,:,Z1_use));
        bwCurr_USE = squeeze(bwNew_vol(:,:,Z1_use));
        if abs(Z1-Z1_use)<round(diskR./2)
            [bwNew, bwCell] = obj.subAddCell(bwCurr_USE, fovImg_USE);
        else
            [bwNew, bwCell] = obj.subAddCell(bwCurr_USE, fovImg_USE,cellMean);
        end
        %% check the plane above/below to make sure no touching..
        Z1_use_tmp = Z1_use + step_Z;
        %first, check if this plane exists
        if (Z1_use_tmp > 0) && (Z1_use_tmp <= nZ)
            bwCell_abovebelow =  squeeze(bwNew_vol(:,:,Z1_use_tmp));
            %check for direct contact between two cells (diagonal
            %through Z doesn't get checked, could dilate to do
            %this)
            tmp = sum(sum(bwCell_abovebelow.*bwCell));
            if tmp > 0
                STOP = 1;
            end
        end
        %%
        if sum(sum(bwCell)) > 0 && STOP == 0
            bwNew_vol(:,:,Z1_use) = bwNew;
            bwNew_vol0(:,:,Z1_use) = bwCell;
            Z1_use = Z1_use + step_Z;
        else
            STOP = 1;
        end
    end
end
%% check and add new spine
index = find(bwNew_vol0);
if ~isempty(index)
    pxNum = length(index);
    obj.saving.changed = 1;
    obj.current.pxNum = pxNum;
    obj.current.X2 = 0;
    obj.current.Y2 = 0;
    obj.cells.nTotal = obj.cells.nTotal+1;
    obj.cells.all(obj.cells.nTotal) = {obj.current};
    obj.cells.labelSpine(index) = obj.cells.nTotal;
    obj.cells.maskSpine(index) = true;
    obj.updateImage();
    fprintf(1, 'Added object #%d: %d pix\n', obj.cells.nTotal, pxNum);
else
    fprintf(1, 'No matching pixels found');
end
    