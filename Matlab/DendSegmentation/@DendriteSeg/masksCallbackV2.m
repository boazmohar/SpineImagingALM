function masksCallback(obj,figH,eventData)
%%
set(obj.handles.ustatH,'String','Choose masks')
pause(0.1)
%% load Sp file for the session
[SpFileName,SpPathName,~] = uigetfile('.mat', 'Select Sp.mat file', 'Sp.mat');
load([SpPathName SpFileName])
%% select FOV stack
[FOVFileName,FOVPathName,~] = uigetfile('.tif', 'Select FOV stack tif');
info = imfinfo([FOVPathName FOVFileName]);
zNum = length(info);
%% setup params
Info = struct();
%%
xSize = info(1).Width;
ySize = info(1).Height;
UM_Zoom = Sp_s.OptionsStruct.UM_1X./Sp_s.OptionsStruct.Zoom;
Pixel_Size = UM_Zoom./Sp_s.OptionsStruct.RefPixels;
imgSize = [xSize, ySize, zNum];

%% For session branches
dend = {};
Ids = Sp_s.All.Id;
uIds = unique(Ids);
uIds(uIds==-1) = [];
%uIds(uIds==1) = [];
[uIdsBranch, ~, ic] = unique(round(uIds));
Session = zeros(xSize, ySize, zNum);
All = zeros(xSize, ySize, zNum);
%% need to take from swc table not All!
for i = 1:length(uIdsBranch)
    idx = find(ic == i);
    branch = round(uIds(idx(1)));
    X1 = Sp_s.Centers{branch}.x;
    Y1 = Sp_s.Centers{branch}.y;
    Z1 = Sp_s.Centers{branch}.z;
    X2 = [];
    Y2 = [];
    Z2 = [];
    for j=1:length(idx)
        currentId = uIds(idx(j));
        subBranch = round((currentId - branch)*10);
        if subBranch == 0
            Ids = find(Sp_s.Table.BranchID == currentId);
            X2 = Sp_s.Table.x(Ids);
            Y2 = Sp_s.Table.y(Ids);
            Z2 = Sp_s.Table.z(Ids);

            % add the privous one if possible
            PrevId = Sp_s.Table.pNum(Ids(1));
            if PrevId ~= -1
                EuclidDistance = sqrt((Sp_s.Table.x(PrevId) - X2(1))^2 + ...
                                      (Sp_s.Table.y(PrevId) - Y2(1))^2 + ...
                                      (Sp_s.Table.z(PrevId) - Z2(1))^2)
                X2 = [Sp_s.Table.x(PrevId); X2];
                Y2 = [Sp_s.Table.y(PrevId); Y2];
                Z2 = [Sp_s.Table.z(PrevId); Z2];
            end
        else
            fields = [0 Sp_s.SubBranch{branch} length(Sp_s.Centers{branch}.x)];
            startField = fields(subBranch);
            endField = fields(subBranch+1);
            X2 = [X2; X1(startField+1:endField)];
            Y2 = [Y2; Y1(startField+1:endField)];
            Z2 = [Z2; Z1(startField+1:endField)];
        end
        
        
    end
    X2 = X2./Pixel_Size;
    Y2 = Y2./Pixel_Size;
    Z2 = Z2./Sp_s.OptionsStruct.Zstep;
    if length(X2) > 1
        
        cumDist=sum(sqrt(diff(X2).^2+diff(Y2).^2+ diff(Z2).^2));
        InterpPoints = round(cumDist)*3;
        V=interparc(InterpPoints,X2,Y2,Z2,'linear');
        [~, index, ~] = unique(round(V),'rows');
        V2 = round(V(sort(index),:));
        V2Index = sub2ind(imgSize,V2(:,1),V2(:,2),V2(:,3));
        Session(V2Index) = 1;
    else
        V2Index = sub2ind(imgSize,round(X2),round(Y2),round(Z2));
        Session(round(V2Index)) = 1;
    end
end
Info.Session = dend;
%%
writetiff(im2uint8(permute(Session, [2, 1, 3])),'Session.tif')
%% trying with actual data from stack
stack = readtiff(FOVPathName, [], FOVFileName);
[sz1, sz2, sz3] = size(stack);
sessionAnatomy = zeros(sz1, sz2, sz3);
%%
X = Sp_s.All.x(Sp_s.TrueFieldsMask==1);
Y = Sp_s.All.y(Sp_s.TrueFieldsMask==1);
Z = Sp_s.All.z(Sp_s.TrueFieldsMask==1);
UM_Zoom = Sp_s.OptionsStruct.UM_1X./Sp_s.OptionsStruct.Zoom;
Pixel_Size = UM_Zoom./Sp_s.OptionsStruct.RefPixels;
X = round(X./Pixel_Size);
Y = round(Y./Pixel_Size);
Z = round(Z./Sp_s.OptionsStruct.Zstep);
pxNum = Sp_s.OptionsStruct.ImagingPixels;
lineNum = Sp_s.OptionsStruct.ImagingLines;
for i =1:length(X)
    sessionAnatomy(Y(i):Y(i)+pxNum-1, X(i):X(i)+lineNum-1, Z(i)) = ...
        stack(Y(i):Y(i)+pxNum-1, X(i):X(i)+lineNum-1, Z(i));
end
%%
writetiff(sessionAnatomy,'SessionAnat.tif')
%% For Session branches
dend = {};
Ids = Sp_s.Table.BranchID;
uIds = unique(Ids);
for i = 1:length(uIds)
    currentId = uIds(i);
    Ids = find(Sp_s.Table.BranchID == currentId);
    X1 = Sp_s.Table.x(Ids);
    Y1 = Sp_s.Table.y(Ids);
    Z1 = Sp_s.Table.z(Ids);
    
    % add the privous one if possible
    PrevId = Sp_s.Table.pNum(Ids(1));
    if PrevId ~= -1
        X1 = [Sp_s.Table.x(PrevId); X1];
        Y1 = [Sp_s.Table.y(PrevId); Y1];
        Z1 = [Sp_s.Table.z(PrevId); Z1];
    end
    X1 = X1./Pixel_Size;
    Y1 = Y1./Pixel_Size;
    Z1 = Z1./Sp_s.OptionsStruct.Zstep;
    X1(X1 > imgSize(1)) = imgSize(1);
    Y1(Y1 > imgSize(2)) = imgSize(2);
    Z1(Z1 > imgSize(3)) = imgSize(3);
    Z1(Z1 < 1) = 1;
    X1(X1 < 1) = 1;
    Y1(Y1 < 1) = 1;
    if length(X1) > 1
        cumDist=sum(sqrt(diff(X1).^2+diff(Y1).^2+ diff(Z1).^2));
        InterpPoints = round(cumDist)*10;
        V=interparc(InterpPoints,X1,Y1,Z1,'linear');
        [~, index, ~] = unique(round(V),'rows');
        V2 = round(V(sort(index),:));
        V2Index = sub2ind(imgSize,V2(:,1),V2(:,2),V2(:,3));
        dend{i} = V2Index;
        All(V2Index) = 1;
    else
        V2Index = sub2ind(imgSize,round(X1),round(Y1),round(Z1));
        dend{i} = V2Index;
        All(V2Index) = 1;
    end
end
Info.All = dend;
%%
writetiff(im2uint8(permute(All, [2,1,3])),'All.tif')
%%
Info.Table = table2array(Sp_s.Table);
[FileName,PathName] = uiputfile
save([PathName FileName], 'Info','-v7')
