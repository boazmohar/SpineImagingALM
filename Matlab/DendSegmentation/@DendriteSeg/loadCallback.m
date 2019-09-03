function loadCallback(obj,figH,eventData)
[filename, pathname] = uigetfile('*.mat', 'Select a mask file');
if isequal(filename,0)
    disp('User selected Cancel')
else
    cells = load(fullfile(pathname, filename));
    %check for correct size
    obj.cells = cells.cells;
    obj.updateImage();
end