%% go to the directory in tier 2
cd('V:\users\Aaron\150814_BMWR17')
load 'Run1ShiftsCorr1'
%% display the correlation of one of the shifts
Shifts = cell2mat( FirstMap(:,1));
Corr = cell2mat( FirstMap(:,2));
X = Shifts(:,1);
Y = Shifts(:,2);
Z = Shifts(:,3 );
close all
figure
scatter3(X,Y,Z,40,Corr,'filled')
colormap hsv
colorbar
[a,b] = max(Corr);ss
title(['x: ' num2str(X(b)) ' y:'  num2str(Y(b)) ' z:'  num2str(Z(b))])

%%
load Run1Try
scatter3(coorList(:,1),coorList(:,2),coorList(:,3),30,Try)
%%
load('Run1newcoorArray')
load('Run1flatTarget')
%%
close all
hFig = figure();
colormap(hFig,'gray');
cameratoolbar();

hAx = axes('Parent',hFig,'Color','black');
view(hAx,37,45);


%
% coorList = newcoorArray;
Zs = unique(coorList(:,3));
for i = 1:length(Zs)
    x = coorList(coorList(:,3)==Zs( i),2);
    x = [min(x) max(x)];
    y = coorList(coorList(:,3)==Zs(i),1);
    y = [min(y) max(y)];
    z = Zs(i);
    [x,y,z] = meshgrid([x(1) x(2)],[y(1) y(2)],z);
   
    hSurf = surface('Parent',hAx,...
        'XData',x,...
        'YData',y,...
        'ZData',z,...
        'FaceColor','texturemap',...
        'CDataMapping','scaled',...
        'FaceLighting','none',...
        'EdgeColor','none');
    data = new(coorList(:,3)==Zs(i));
      set(hSurf,'CData',data)
      pause
end
%%
scatter3(coorList(:,1),coorList(:,2),coorList(:,3),30,flatTarget)
%%
flatTarget2= reshape(Mean,105840,1);
%%
Index1 = 1:72;
Index2 = 1:35;
[X,Y] = meshgrid(Index1,Index2);
X = X(:);
Y= Y(:);
%%
figure
scatter3(X,Y,coorList(1:2520,3),30,flatTarget2(1:2520))
figure()
imshow(Mean(:,:,1),[])
%% unpacking the mean from flat 1d array
Zs = unique(coorList(:,3));
Volume = zeros(80,200,length(Zs));
Counter = zeros(80,200,length(Zs));
for i = 1: length(Zs)
    Z =Zs(i); 
    Zindexs = Z == coorList(:,3);
    X = coorList(Zindexs,1);
    Y = coorList(Zindexs,2);
    Data = newMean2(Zindexs);
    X = X - min(X) +1;
    Y = Y - min(Y) + 1;
    for j = 1:length(X)
        Volume(X(j),Y(j),i) = Volume(X(j),Y(j),i)*Counter(X(j),Y(j),i) + ...
            Data(j);
        Counter(X(j),Y(j),i) = Counter(X(j),Y(j),i)+1;
    end
end
%
writetiff(Volume,'newRegMean2.tif')