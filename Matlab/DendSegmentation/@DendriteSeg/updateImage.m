function updateImage(obj)
Z = obj.display.Z;
plane = obj.display.volume(:,:,Z);
scaledImage = obj.normImage(plane);
mask = obj.cells.maskSpine(:,:,Z);
mask2 = obj.cells.maskBorder(:,:,Z);
if obj.display.showDend
    mask3 = obj.cells.maskDend(:,:,Z);
    sahdedImage = imShade(scaledImage,mask,mask2,mask3);
else
    sahdedImage = imShade(scaledImage,mask,mask2);
end
set(obj.handles.imH,'CData',sahdedImage)
%% dendrite numbering
if isfield(obj.display,'dendriteMedian')
    dendZstart = obj.display.dendriteMedian(:,3);
    dendZend = obj.display.dendriteMedian(:,6);
    for i = 1:length(dendZstart)
        if Z > dendZstart(i) - 2 && Z < dendZend(i) + 2
            set(obj.display.dendText(i,1),'visible','on')
            set(obj.display.dendText(i,2),'visible','on')
        else
            set(obj.display.dendText(i,1),'visible','off')
            set(obj.display.dendText(i,2),'visible','off')
        end
    end
end
    
