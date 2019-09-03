function Error  = getMeta(obj)
%% disk
Error = 0;
 diskR = str2double(get(obj.handles.udH, 'String'));
if isnan(diskR)
    fprintf(1, 'Error reading disk radius: %s, try again\n', ...
        get(obj.handles.udH, 'String'));
    Error = 1;
    return
end
if (diskR < 1) || diskR > 50
    fprintf(1, 'Disk too small or too big, try again\n');
    Error = 1;
    return
end
if ~iswithintol(round(diskR), diskR, 10^2*eps)
    fprintf(1, 'Disk radius must be an integer, try again\n');
    Error = 1;
    return
end

%% threshold
cThresh = str2double(get(obj.handles.utH, 'String'));
if isnan(cThresh)
    fprintf(1, 'Error reading threshold: %s, try again\n', ...
        get(obj.handles.utH, 'String'));
    Error = 1;
    return
end
if cThresh <= 0 || cThresh >= 100
    fprintf(1, 'Threshold too small or too big, try again\n');
    Error = 1;
    return
end
%% dendrite number

dendNum = str2double(get(obj.handles.unH, 'String'));
if isnan(dendNum)
    fprintf(1, 'Error reading dendNum: %s, try again\n', ...
        get(obj.handles.unH, 'String'));
    Error = 1;
    return
end
if dendNum <= 0 || dendNum >= 100
    fprintf(1, 'dendNum too small or too big, try again\n');
end

if ~iswithintol(round(dendNum), dendNum, 10^2*eps)
    fprintf(1, 'dendNum must be an integer, try again\n');
    Error = 1;
    return
end
%% Quality
Value = get(obj.handles.uqH, 'Value');
Quality = get(obj.handles.uqH,'String');
Quality = Quality{Value};
%% Note
Note = get(obj.handles.unoteH,'String');
%% building current selection
obj.current.diskR = diskR;
obj.current.cThresh = cThresh ;
obj.current.Quality = Quality;
obj.current.dendNum = dendNum;
obj.current.note = Note;
%% se2
obj.display.se2 = strel('disk',round(diskR),4);