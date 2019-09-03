function addBorder(obj)
set(obj.handles.ustatH,'String','Add border')
%% get a second point
X = obj.current.X;
Y = obj.current.Y;
Z = obj.current.Z;
[X2,Y2] = ginputc(1,'color',[1 1 1],'LineWidth',1);
X2 = round(X2);
Y2 = round(Y2);
if obj.cells.maskSpine(Y2,X2,Z)
    disp('** Trying to add a border in existing cell region, try again');
    return
end
%% make a line in XY
allX = [X X2];
allY = [Y Y2];
cumDist=sum(sqrt(diff(allX).^2+diff( allY).^2));
InterpPoints = round(cumDist^2);
V=interparc(InterpPoints,[X X2],[Y Y2],'linear');
V2 = unique(round(V),'rows');
V2(:,3) =repmat(Z,size(V2,1),1);
%% expand it in Z +-5 planes
Len = size(V2,1);
V3 = repmat(V2,19,1);
Vz = meshgrid(-9:9,1:Len);
Vz = Vz(:);
V3(:,3) = V3(:,3) + Vz;
badLines = V3(:,3) < 1 | V3(:,3) >= obj.display.nZ;
V3(badLines,:) = [];
V3index = sub2ind(size(obj.cells.maskSpine),V3(:,2),V3(:,1),V3(:,3));
% tempVol = false(size(obj.cells.maskSpine));
% tempVol(V3index) = true;
%% check if there is a mask  there an add
if sum(obj.cells.maskSpine(V3index)) > 0
    disp('Border is touching a allready masked area')
    return;
else
    obj.saving.changed = 1;
    obj.cells.nTotal = obj.cells.nTotal+1;
    obj.current.X2 = X2;
    obj.current.Y2 = Y2;
    obj.cells.labelBorder(V3index) = obj.cells.nTotal;
    obj.cells.maskBorder(V3index) = true;
    obj.current.pxNum = length(V3index);
    obj.cells.all(obj.cells.nTotal) = {obj.current};
    obj.cells.labelSpine(V3index) = obj.cells.nTotal; 
    obj.cells.maskSpine(V3index) = true;
    obj.updateImage();
    fprintf(1, 'Added border #%d: %d pix\n', obj.cells.nTotal, obj.current.pxNum); 
end