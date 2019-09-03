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
lastFOV = [];
%% get all expanded*.tif full paths
for d=1:length(P)
    p1 = P{d};
    cd(p1)
    disp(p1);
    f1 = F{d};
    fSp = cellfun(@(x) ~isempty(strfind(x, 'Sp')), f1);
    fSp2 = cellfun(@(x) ~isempty(strfind(x, 'Sp_2')), f1);
    if sum(fSp2)
        disp('skipped')
        continue;
    end
    if sum(fSp) >= 1
        disp('found')
        %% load
        a = load('Sp');
        Sp_s = a.Sp_s;
        Sp_s.tableArray = table2array(Sp_s.Table);
        save('Sp_2', 'Sp_s')
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
