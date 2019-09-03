function BinaryMaskFromSp(obj)
%% load Sp file for the session
[SpFileName,SpPathName,~] = uigetfile('.mat', 'Select Sp.mat file', 'Sp.mat');
load([SpPathName SpFileName])
%% select FOV stack
[FOVFileName,FOVPathName,~] = uigetfile('.tif', 'Select FOV stack tif');
info = imfinfo([FOVPathName FOVFileName]);
zNum = length(info);
FOV = readtiff([FOVPathName FOVFileName]);
%% select Reconstraction Stack
[ReconFileName,ReconPathName,~] = uigetfile('.tif', 'Select Reconstraction stack tif');
Stack = readtiff([ReconPathName ReconFileName]);
%% setup params
Info = struct();
Info.FOV = FOV;
Info.Stack = Stack;

UM_Zoom = Sp_s.OptionsStruct.UM_1X./Sp_s.OptionsStruct.Zoom;
Pixel_Size = UM_Zoom./Sp_s.OptionsStruct.RefPixels;
imgSize = [Sp_s.OptionsStruct.RefPixels, Sp_s.OptionsStruct.RefPixels, zNum];

%% For session branches
dend = {};
Ids = Sp_s.All.Id;
uIds = unique(Ids);
uIds(uIds==-1) = [];
uIds(uIds==1) = [];
for i = 1:length(uIds)
    currentId = uIds(i);
    X1 = Sp_s.All.x(Ids == currentId)./Pixel_Size;
    Y1 = Sp_s.All.y(Ids == currentId)./Pixel_Size;
    Z1 = Sp_s.All.z(Ids == currentId)./Sp_s.OptionsStruct.Zstep;
    if length(X1) > 1
        
        cumDist=sum(sqrt(diff(X1).^2+diff(Y1).^2+ diff(Z1).^2));
        InterpPoints = round(cumDist)*3;
        V=interparc(InterpPoints,X1,Y1,Z1,'linear');
        [~, index, ~] = unique(round(V),'rows');
        V2 = round(V(sort(index),:));
        V2Index = sub2ind(imgSize,V2(:,1),V2(:,2),V2(:,3));
        dend(i) = {V2Index};
    else
        V2Index = sub2ind(imgSize,X1,Y1,Z1);
        dend(i) = {V2Index};
    end
end
Info.Session = dend;
%% For Session branches
dend = {};
Ids = Sp_s.Table.BranchID;
uIds = unique(Ids);
for i = 1:length(uIds)
    currentId = uIds(i);
    X1 = Sp_s.Table.x(Ids == currentId)./Pixel_Size;
    Y1 = Sp_s.Table.y(Ids == currentId)./Pixel_Size;
    Z1 = Sp_s.Table.z(Ids == currentId)./Sp_s.OptionsStruct.Zstep;
    if length(X1) > 1
        cumDist=sum(sqrt(diff(X1).^2+diff(Y1).^2+ diff(Z1).^2));
        InterpPoints = round(cumDist)*10;
        V=interparc(InterpPoints,X1,Y1,Z1,'linear');
        [~, index, ~] = unique(round(V),'rows');
        V2 = round(V(sort(index),:));
        V2Index = sub2ind(imgSize,V2(:,1),V2(:,2),V2(:,3));
        dend(i) = {V2Index};
    else
        V2Index = sub2ind(imgSize,round(X1),round(Y1),round(Z1));
        dend(i) = {V2Index};
    end
end
Info.All = dend;
%%
Info.Table = Table2array(Sp_s.Table);
[FileName,PathName] = uiputfile
save([PathName FileName], 'Info','-v7')
