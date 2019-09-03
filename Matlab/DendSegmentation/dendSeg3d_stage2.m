%Transform registered mean of mROI data to sparse interpolated volume

pixelSize=[0.4 0.4 1];
PSF = [.4 .4 2.5];
expdir='V:\users\Aaron\Database\BMWR30\FOV1\151111Run1';
expname='expended';





%%
embeddedFields=readtiff(expdir,[],[expname,'_new']);
%%
cell_vol = embeddedFields;
cell_vol(isnan(cell_vol))=0;
si=size(cell_vol);
%upsample
A=eye(4)*2;
A(4,4)=1;
R=makeresampler('cubic','bound');
tform=maketform('affine',A);
tic
cell_vol_up = tformarray(cell_vol, tform, R, [1 2 3], [1 2 3], [si(1)*2 si(2)*2, si(3)*2],[],0);
toc
%
tic
filt_vol_up = imgaussfilt3(cell_vol_up,PSF./pixelSize/2.355*2);
toc
top_vol_up = filt_vol_up;
%
% pixelSize=pixelSize./2;
% tic
% top_vol_up=imtophat(filt_vol_up,distStrel3D(0.5, pixelSize));
% toc        
%
% labelimg_up = imCellEditInteractive_3DMA_110407MA(cell_vol_up,[],[],1,1);
%%
cd(expdir)
save('top_vol_up','top_vol_up')
%%
cd(expdir)
load('top_vol_up')
%%

dend = DendriteSeg(top_vol_up);
%%

all = dend.cells.all;
for i=1:length(all)
    if (all{i}.dendNum == 4)
        dend.cells.all{i}.dendNum = 2;
        disp(i)
    end
end
%%
writetiff(top_vol_up,'top_vol_up')
%%
RegInfo=dendSeg3d(RegInfo,embeddedFields, pixelSize);

save('X:\users\Aaron\150821_BMWR17\Run1labelimg.mat','labelimg');

%
load('d:\data\BMWR30\mask2.mat')
%
load('G:\Scripts\GeneralMat\DendSegmentation\mask3.mat')
%
embeddedFields2=permute(embeddedFields,[2 1 3]);
Remove spines and replace values with nearby surround 
labelimg=permute(cells.downsampled.labelimg,[2 1 3]);
spineRemDend=embeddedFields.*~(imdilate(labelimg>0,distNhood3D(0.5, pixelSize)));
spineRemDend(isnan(spineRemDend))=0;
spineRemDend2=spineRemDend;
for i=1:max(labelimg(:))
    surroundMask=imdilate(labelimg==i,distNhood3D(2, pixelSize))-imdilate(labelimg==i,distNhood3D(0.5, pixelSize));
    replaceVal=nanmean(spineRemDend(surroundMask>0));
    spineRemDend2(labelimg==i)=replaceVal;
end    
%
smooth removed spine regions
smoothDend=imgaussfilt3(spineRemDend2,0.75);
%
prerform tophat with SE about dendrite size
smoothDendTH=imtophat(smoothDend,distStrel3D(2, pixelSize));


writetiff(smoothDend,'G:\Scripts\GeneralMat\DendSegmentation\smoothDend.tif')
writetiff(smoothDendTH,'G:\Scripts\GeneralMat\DendSegmentation\smoothDendTH.tif')

%
%Manually fix incorrectly merged dendrite by cleaving in ImageJ

%%
smoothDendTH=readtiff('V:\users\Aaron\150821_BMWR17\',[],'smoothDendTH.tif');
testvol=smoothDendTH>Thresh;
%skeletonize and display
Thresh=16;
skel = Skeleton3D(smoothDendTH>Thresh);

w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% convert skeleton to graph structure
[A,node,link] = Skel2Graph3D(skel,10);

% convert graph structure back to (cleaned) skeleton
skel = Graph2Skel3D(node,link,w,l,h);

writetiff(double(skel),'X:\users\Aaron\150821_BMWR17\dendSkel.tif')

figure();
col=[.7 .7 .8];
hiso = patch(isosurface(testvol,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(testvol,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(testvol,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;
w=size(skel,1);
l=size(skel,2);
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
plot3(y,x,z,'square','Markersize',4,'MarkerFaceColor','r','Color','r');            
set(gcf,'Color','white');
view(140,80)
%%

w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% convert skeleton to graph structure
[A,node,link] = Skel2Graph3D(skel,1);

% convert graph structure back to (cleaned) skeleton
skel2 = Graph2Skel3D(node,link,w,l,h);



% iteratively convert until there are no more 2-nodes left
[A2,node2,link2] = Skel2Graph3D(skel2,10);
% while(min(cellfun('length',{node2.conn}))<3)
%     skel2 = Graph2Skel3D(node2,link2,w,l,h);
%     [A2,node2,link2] = Skel2Graph3D(skel2,4);
% end;

% display result
figure();
hold on;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    for j=1:length(node2(i).links)    % draw all connections of each node
        if(node2(i).conn(j)<1)
            col='b'; % branches are blue
        else
            col='r'; % links are red
        end;
        
        % draw edges as lines using voxel positions
        for k=1:length(link2(node2(i).links(j)).point)-1            
            [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',3);
        end;
    end;
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor','y',...
        'Color','k');
end;
axis image;axis off;
set(gcf,'Color','white');
drawnow;
view(-17,46);

%%
%Create dendrite labelimg
load('V:\users\Aaron\150821_BMWR17\Run1labelimg.mat');
% labelimg=RegInfo.Cell.labelimg;
um2bin=10;

dendLabelImg=zeros(size(labelimg));
dendLabelTable=zeros(1,2);
counter=max(labelimg(:))+1;
for i=1:length(link2)
    [X,Y,Z]=ind2sub(size(labelimg),link2(i).point);
    dist=sqrt((diff(X).*pixelSize(1)).^2+(diff(Y).*pixelSize(2)).^2+(diff(Z).*pixelSize(3)).^2);
    tDist=sum(dist);
    nSeg=round(tDist./um2bin);
    stopPt=0;
    for j=1:nSeg
        dendLabelImgTemp=zeros(size(labelimg));
        startPt=stopPt+1;
        stopPt=min([startPt+find(cumsum(dist(startPt:end))>um2bin,1,'first'),length(link2(i).point)]);
        dendLabelImgTemp(link2(i).point(startPt:stopPt))=1;
        dendLabelImgTemp=imdilate(dendLabelImgTemp>0,distNhood3D(1, pixelSize));
        dendLabelImg(dendLabelImgTemp)=counter;
        dendLabelTable(counter,1)=i;
        dendLabelTable(counter,2)=j;
        counter=counter+1;  
    end   
end



labelimg2=labelimg;
for i=1:max(labelimg(:))
    labelimgTemp=zeros(size(labelimg));
    labelimgTemp(labelimg==i)=1;
    labelimgTemp=imdilate(labelimgTemp,distNhood3D(0.5, pixelSize));
    labelimg2(labelimgTemp==1)=i;    
end    

dendLabelImg=dendLabelImg.*~(imdilate(labelimg2>0,distNhood3D(0.5, pixelSize)));
labelimgAll=dendLabelImg+labelimg;
save('X:\users\Aaron\150821_BMWR17\Run1labelimgAll.mat','labelimgAll','-v7.3');

writetiff(labelimgAll,'X:\users\Aaron\150821_BMWR17\Run1labelimgAll.tif')

%%
load([expdir,expname,'labelimgAll.mat'])
labelimgAll=permute(labelimgAll,[2 1 3]);
linearLabelimg=labelimgAll(:);
nLabel=size(dendLabelTable,1);
labelInfoAll=cell(1,length(fieldsTform));
for i=1:length(fieldsTform)
    Tform=fieldsTform{i};
    labelInfo=zeros(size(Tform,1),3);
    for j=1:size(Tform,1)
        labelInfo(j,1)=linearLabelimg(Tform(j,1));
        labelInfo(j,2)=Tform(j,2);
        labelInfo(j,3)=i;
    end
    labelInfoAll{i}=labelInfo';
end 
labelArray=cell2mat(labelInfoAll);

labelCell=cell(1,nLabel);
for i=1:nLabel
    labelCell{i}=labelArray(2:3,labelArray(1,:)==i)';
end
save([expdir,expname,'labelCell.mat'],'labelCell')

%%
branchOrder=[10,9,6,1,2];
dir=[-1,-1,-1,-1,1];
orderedDendIdx=[];
for i=1:length(branchOrder)
    if dir(i)==1
        orderedDendIdx=[orderedDendIdx,find(dendLabelTable(:,1)==branchOrder(i))'];
    else
        orderedDendIdx=[orderedDendIdx,flip(find(dendLabelTable(:,1)==branchOrder(i)))'];
    end    
end    

badSegments=[1,13]; %not in actual mROI space
orderedDendIdx(badSegments)=[];

%%
branchOrder=[8];
dir=[1];
orderedDendIdx1p1=[];
for i=1:length(branchOrder)
    if dir(i)==1
        orderedDendIdx1p1=[orderedDendIdx1p1,find(dendLabelTable(:,1)==branchOrder(i))'];
    else
        orderedDendIdx1p1=[orderedDendIdx1p1,flip(find(dendLabelTable(:,1)==branchOrder(i)))'];
    end    
end    

badSegments=[7]; %not in actual mROI space
orderedDendIdx1p1(badSegments)=[];
%%
dir=[-1,1];
branchOrder2=[3,4];
orderedDendIdx2=[];
for i=1:length(branchOrder2)
    orderedDendIdx2=[orderedDendIdx2,find(dendLabelTable(:,1)==branchOrder2(i))'];
end    
badSegments=[3,10]; %not in actual mROI space, crossing point
orderedDendIdx2(badSegments)=[];
