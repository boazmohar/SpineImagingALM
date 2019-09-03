function scrollCallback(obj, ~, eventData)

value = eventData.VerticalScrollCount;
currentZ = obj.display.Z;
nextZ = currentZ + value;
if nextZ > 0 && nextZ < obj.display.nZ
    set(obj.handles.uzH,'string',num2str(nextZ))
    obj.display.Z = nextZ;
    obj.updateImage()
end