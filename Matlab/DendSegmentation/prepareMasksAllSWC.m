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
    fAuto = cellfun(@(x) ~isempty(strfind(x, 'prepareMasksAutoSWC')), f1);
    if sum(fAuto)
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
 
        UM_Zoom = Sp_s.OptionsStruct.UM_1X./Sp_s.OptionsStruct.Zoom;
        Pixel_Size = UM_Zoom./Sp_s.OptionsStruct.RefPixels;
        imgSize = [sz1, sz2, sz3];
     
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
        Info.AllSWC = dend;
        Info.TableSWC = table2array(Sp_s.Table);
        %         AllInfo(d) = {Info};
        save('prepareMasksAutoSWC.mat', 'Info','-v7');
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
