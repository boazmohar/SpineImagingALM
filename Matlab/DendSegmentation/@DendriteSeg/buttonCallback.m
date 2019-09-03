    function buttonCallback(obj,figH,~)
% button click callback handles adding and deleting spines

%% get current location
Position =obj.handles.pxInfoH;
PosText = get(Position, 'String');
C = strsplit(PosText,{'(',')',',' ' '});
X = str2num(C{2});
Y = str2num(C{3});
%% error checking and getting current meta data
if isempty(X) || isempty(Y)
    return;
end
Error = obj.checkLimits(X,Y);
Error = Error || obj.getMeta();
obj.current.X = X;
obj.current.Y = Y;
obj.current.Z = obj.display.Z;
%% sdd or remove mask
Z = obj.display.Z;
if ~Error
    Type = get(figH,'SelectionType');
    if obj.cells.maskSpine(Y,X,Z) == 1 % if clicked location has a mask 
        if strcmpi(Type,'extend')
            % delete current as the click was with Shift
            obj.deleteMask(X,Y,Z);
        else
            % normal click on mask location -> error
            disp('** Trying to add in existing cell region, try again');
        end
    else
        % no mask in the current location
        if strcmpi(get(figH,'SelectionType'),'extend')
            % if Ctrl click --> error
            fprintf(1, '** Trying to delete, not on a cell, try again\n');
        elseif strcmpi(Type,'alt')
            % add a border as this was a Ctrl+click
            obj.addBorder();
        else
            % try to add a new spine mask
            obj.addCell3D()
        end
        
    end
else
    disp('Error adding cell')
end

set(obj.handles.ustatH,'String','Ready')


