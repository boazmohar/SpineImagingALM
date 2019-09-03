 %% upsample and top hat all the expanded images
% close all; clear; clc;
%% prepare
pixelSize=[0.4 0.4 1];
PSF = [.4 .4 2.5];
pixelSize2=pixelSize./2;
A=eye(4)*2;
A(4,4)=1;
R=makeresampler('cubic','bound');
tform=maketform('affine',A);
structure = distStrel3D(0.5, pixelSize2);
%% all directories under database
[a, name] = system('hostname');
if strfind(name, 'moharb-lm1')
    path = '/Volumes/svobodalab/Database';
else
    path = 'V:\users\Aaron\Database';
end
[P,F] = subdir(path);
% all directories that have Run
RunIndex = cellfun(@(x) ~isempty(strfind(lower(x), 'run')) && ...
    isempty(strfind(lower(x), 'old')),P);
P = P(RunIndex);
F = F(RunIndex);
%% get all expanded*.tif full paths
for i=1:length(P)
    p1 = P{i};
    cd(p1)
    
    disp(p1);
    f1 = F{i};
    fIndex = cellfun(@(x) ~isempty(strfind(x, 'exp')), f1);
    ftop = cellfun(@(x) ~isempty(strfind(x, 'top')), f1);
    fup2 = cellfun(@(x) ~isempty(strfind(x, 'up2')), f1);
    ffilt = cellfun(@(x) ~isempty(strfind(x, 'filt_vol_up')), f1);
    if ~sum(fIndex) || (sum(ftop) && (sum(fup2) || sum(ffilt))) 
        disp('skipped')
        continue;
    end
    f = f1{fIndex};  
    cell_vol=readtiff(p1,[],f);
    cell_vol(isnan(cell_vol))=0;
    si=size(cell_vol);
    disp('upsample')
    tic
    cell_vol_up = tformarray(cell_vol, tform, R, [1 2 3], [1 2 3], ...
                             [si(1)*2 si(2)*2, si(3)*2],[],0);
    toc
    %
    disp('filter')
    tic
    filt_vol_up = imgaussfilt3(cell_vol_up,PSF./pixelSize/2.355*2);
    toc
    save('filt_vol_up','filt_vol_up', '-v7.3')
    writetiff(filt_vol_up,'filt_vol_up')
    
    disp('tophat')
    tic
    top_vol_up=imtophat(filt_vol_up,structure);
    toc 
    
    save('top_vol_up','top_vol_up', '-v7.3')
    writetiff(top_vol_up,'top_vol_up')
    disp('saved')
end
