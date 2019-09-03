
%% upsample and top hat all the expanded images
% close all; clear; clc;

%% all directories under database
path = 'V:\users\Aaron\Database';
[P,F] = subdir(path);
% all directories that have Run
RunIndex = cellfun(@(x) ~isempty(strfind(lower(x), 'run')) && ...
    isempty(strfind(lower(x), 'old')),P);
P = P(RunIndex);
F = F(RunIndex);
AllInfo = cell(1,length(P));
%% get all expanded*.tif full paths
for d=1:length(P)
    p1 = P{d};
    cd(p1)
    disp(p1);
    f1 = F{d};
    fIndex = cellfun(@(x) ~isempty(strfind(x, 'exp')), f1);
    fSp = cellfun(@(x) ~isempty(strfind(x, 'Sp')), f1);
    fAuto = cellfun(@(x) ~isempty(strfind(x, 'prepareMasksAuto')), f1);
    fAuto2 = cellfun(@(x) ~isempty(strfind(x, 'prepareMasksAutoSWC')), f1);
    if sum(fAuto & not(fAuto2))
        disp('skipped')
        continue;
    end
    if sum(fIndex) + sum(fSp) >= 2
        disp('found')
        %% load
        a = load('Sp');
        Sp_s = a.Sp_s;
        Ids = Sp_s.All.Id;
        current = pwd;
        cd('..')
        stack = readtiff(pwd, [], 'Stack_32bit.tif');
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
        %uIds(uIds==1) = [];
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
                    ZIndex = round(V2(k,3));
                    ZIndex(ZIndex > sz3) = sz3;
                    ZIndex(ZIndex < 1) = 1;
                    Session(YIndex, XIndex, ZIndex) = stack(YIndex, XIndex, ZIndex);
                end
            else
                YIndex = round(Y2)-before1:round(Y2)+after1;
                XIndex = round(X2)-before2:round(X2)+after2;
                
                YIndex(YIndex > sz1) = [];
                YIndex(YIndex < 1) = [];
                XIndex(XIndex > sz2) = [];
                XIndex(XIndex < 1) = [];
                ZIndex = round(Z2);
                ZIndex(ZIndex > sz3) = sz3;
                ZIndex(ZIndex < 1) = 1;
                Session(YIndex, XIndex, ZIndex) = stack(YIndex, XIndex, ZIndex);
                
            end
        end
        %%
        writetiff(imgaussfilt3(Session,1),'Session.tif')
        %%
        dend = {};
        Ids = Sp_s.Table.BranchID;
        uIds = unique(Ids);
        for i = 1:length(uIds)
            currentId = uIds(i);
            Ids = find(Sp_s.Table.BranchID == currentId);
            X1 = Sp_s.Table.x(Ids);
            Y1 = Sp_s.Table.y(Ids);
            Z1 = Sp_s.Table.z(Ids);
            
            %add the privous one if possible
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
            else
                V2Index = sub2ind(imgSize,round(X1),round(Y1),round(Z1));
                dend{i} = V2Index;
            end
        end
        Info.All = dend;
        Info.Table = table2array(Sp_s.Table);
        %% redo with swc file
        [pathstr,name,ext]= fileparts( Sp_s.OptionsStruct.filename);
        current = pwd;
        cd('..')
        swcFile = [pwd '\' name ext];
        cd(current)
        [fid, ~] = fopen(swcFile);
        if fid==-1
            throw(MException('SpineImaging:readSWC','Error opening file %s'))
        end
        tline = fgets(fid);
        LineIgnore =0;
        while ischar(tline) && tline(1) == '#'
            LineIgnore = LineIgnore+1;
            tline = fgets(fid);
        end
        fclose(fid);
        %% read the table
        VarNames = {'Num','Type','x','y','z','r','pNum'};
        Table = readtable(swcFile,'FileType','text','ReadVariableNames',0,'Delimiter',' ',...
            'HeaderLines',LineIgnore);
        Table.Properties.VariableNames = VarNames;
        %% Converting pixels to um
        UM_Zoom = Sp_s.OptionsStruct.UM_1X./Sp_s.OptionsStruct.Zoom;
        Pixel_Size = UM_Zoom./Sp_s.OptionsStruct.RefPixels;
        Table.x = Table.x.*Pixel_Size;
        Table.y = Table.y.*Pixel_Size;
        Table.z = Table.z.*Sp_s.OptionsStruct.Zstep;
        %% setting the object's table
        Sp_s.Table = Table;
        %% getting branch number
        Base = find(Sp_s.Table.pNum == -1);
        % if isempty(Base) || length(Base) > 1
        %     throw(MException('SpineImaging:c_DetectBranches','Error with the tree: check disconnected branches'));
        % end
        %% init branch id coulmn of the table
        Branchs = zeros(height(Sp_s.Table),1);
        Sp_s.Table.BranchID = Branchs;
        BranchID = 0;
        %% go over the tree
        while ~isempty(Base)
            BranchID = BranchID+1;
            Sp_s.Table.BranchID(Base) = BranchID;
            Parents = find(Sp_s.Table.pNum == Base(1));
            if isempty(Parents)
                disp(['initPath:base - Error with the tree: Only one node: #' num2str(Base(1))]);
                Base(1) = [];
                continue;
            elseif length(Parents) > 1
                Base(1) = [];
                Base = [Base; Parents(2:end)];
                BranchID = BranchID+1;
            else
                Base(1) = [];
            end
            %% go over the rest of the branch
            Current = Parents(1);
            Sp_s.Table.BranchID(Current) = BranchID;
            Flag = 1;
            while Flag
                Parents = find(Sp_s.Table.pNum == Current);
                if isempty(Parents)
                    Flag = 0;
                elseif length(Parents) > 1
                    Base = [Base; Parents(2:end)];
                    Current = Parents(1);
                    BranchID = BranchID+1;
                    Sp_s.Table.BranchID(Current) = BranchID;
                else
                    Current = Parents(1);
                    Sp_s.Table.BranchID(Current) = BranchID;
                end
            end
        end
        Sp_s.BranchNum = BranchID;
        %%
        %         AllInfo(d) = {Info};
        save('prepareMasksAuto.mat', 'Info','-v7');
        %
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
