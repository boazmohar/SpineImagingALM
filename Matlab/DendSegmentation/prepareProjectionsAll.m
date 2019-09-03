%% upsample and top hat all the expanded images
close all; clear; clc;

%% all directories under database
path = 'V:\users\Aaron\Database';
[P,F] = subdir(path);
% all directories that have Run
RunIndex = cellfun(@(x) ~isempty(strfind(lower(x), 'run')) && ...
    isempty(strfind(lower(x), 'old')),P);
P = P(RunIndex);
F = F(RunIndex);
AllInfo = cell(1,length(P));
lastFOV = [];
%% get all expanded*.tif full paths
for d=1:length(P)
    p1 = P{d};
    cd(p1)
    disp(p1);
    f1 = F{d};
    fIndex = cellfun(@(x) ~isempty(strfind(x, 'exp')), f1);
    fSp = cellfun(@(x) ~isempty(strfind(x, 'Sp')), f1);
    fAuto = cellfun(@(x) ~isempty(strfind(x, 'project')), f1);
%     if sum(fAuto)
%         disp('skipped')
%         continue;
%     end
    if sum(fIndex) + sum(fSp) >= 2
        disp('found')
        %% load
        a = load('Sp');
        Sp_s = a.Sp_s;
        Ids = Sp_s.All.Id;
        current = pwd;
        cd('..')
        if(isempty(lastFOV) || ~strcmpi(pwd, lastFOV))
            lastFOV = pwd;
            stack = readtiff(pwd, [], 'Stack_32bit.tif');
        end
        cd(current)
        %% setup params
        Info = struct();
        [sz1, sz2, sz3] = size(stack);
%         zNum = sz3;
%         xSize = sz2;
%         ySize = sz1;
        UM_Zoom = Sp_s.OptionsStruct.UM_1X./Sp_s.OptionsStruct.Zoom;
        Pixel_Size = UM_Zoom./Sp_s.OptionsStruct.RefPixels;
        imgSize = [sz1, sz2, sz3];
%         before1 = 14; %27
%         after1 = 13; %26
%         before2 = 27; %54
%         after2 = 26; % 53
        before1 = 27;
        after1 = 26;
        before2 = 54;
        after2 = 53;
        %% For session branches
        Ids = Sp_s.All.Id;
        uIds = unique(Ids);
        uIds(uIds==-1) = [];
        uIds = [1 ;uIds];
        [uIdsBranch, ~, ic] = unique(floor(uIds));
        Session = zeros(sz1, sz2, sz3);
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
                            (Sp_s.Table.z(PrevId) - Z2(1))^2);
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
                for k=1:size(V2,1)
                    YIndex = round(V2(k,2))-before1:round(V2(k,2))+after1;
                    XIndex = round(V2(k,1))-before2:round(V2(k,1))+after2;
                    YIndex(YIndex > sz1) = [];
                    YIndex(YIndex < 1) = [];
                    XIndex(XIndex > sz2) = [];
                    XIndex(XIndex < 1) = [];
                    Session(YIndex, XIndex, round(V2(k,3))) = stack(YIndex, XIndex, round(V2(k,3)));
                end
            else
                YIndex = round(Y2)-before1:round(Y2)+after1;
                XIndex = round(X2)-before2:round(X2)+after2;
                
                YIndex(YIndex > sz1) = [];
                YIndex(YIndex < 1) = [];
                XIndex(XIndex > sz2) = [];
                XIndex(XIndex < 1) = [];
                if branch==1
                    cell =  stack(YIndex, XIndex, round(Z2));
                    cell = cell - min(cell(:));
                    cell = (cell / max(cell(:))) * 10;
                    Session(YIndex, XIndex, round(Z2)) = cell;
                else
                    Session(YIndex, XIndex, round(Z2)) = stack(YIndex, XIndex, round(Z2));
                end
                
            end
        end
        %%
        XY = nanmean(Session,  3);
        XY = XY - min(XY(:));
        XY = XY / max(XY(:));
        XZ = squeeze(nanmean(Session,  1));
        XZ = XZ - min(XZ(:));
        XZ = XZ / max(XZ(:));
        YZ = squeeze(nanmean(Session,  2));
        YZ = YZ - min(YZ(:));
        YZ = YZ / max(YZ(:));
        All1 = [XY YZ];
        square = nan(size(XZ,2));
        XZ2 = [XZ; square];
        All1 = [All1; XZ2'];
        imwrite(All1, 'projections2.png')
        filename = current(25:end);
        filename = strrep(filename, '\','_');
        path = current(1:24);
        full = strcat(path ,'projections\', filename ,'_projections.png');
        imwrite(All1, full)
        %%
    else
        disp('skipped')
        continue;
    end
    
end

%%
% for i=1: length(AllInfo)
%     if ~isempty(AllInfo(i))
%         p1 = P{i};
%         cd(p1)
%         disp(p1);
%         Info = AllInfo(i);
%         save('prepareMasksAuto.mat', 'Info','-v7');
%     end
% end
