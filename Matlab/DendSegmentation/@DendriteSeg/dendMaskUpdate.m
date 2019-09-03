function dendMaskUpdate(obj)
set(obj.handles.ustatH,'String','Mask Update')
Sz = size(obj.cells.maskSpine);
obj.cells.maskDend = false(Sz);
obj.display.dendriteMedian = zeros(obj.cells.BranchNum,6);
obj.cells.branches = {};
for i = 1:obj.cells.BranchNum
    Ids = find(obj.cells.dendriteTable.BranchID == i);
    if length(Ids) <= 2
        X1 = obj.cells.dendriteTable.y(Ids(1));
        Y1 = obj.cells.dendriteTable.x(Ids(1));
        Z1 = obj.cells.dendriteTable.z(Ids(1));
        obj.display.dendriteMedian(i,:) = [X1, Y1, Z1 ,X1, Y1, Z1 ];
        obj.display.dendText(i,1) = text(Y1-10, X1-10,num2str(i), ...
            'fontsize',14,'color',[1 1 1],'visible','off');
        obj.display.dendText(i,2) = text(Y1+10, X1+10,num2str(i), ...
            'fontsize',14,'color',[1 1 1],'visible','off');
        V2 = [Y1, X1, Z1];
        if length(Ids) == 2
            X1 = obj.cells.dendriteTable.y(Ids(2));
            Y1 = obj.cells.dendriteTable.x(Ids(2));
            Z1 = obj.cells.dendriteTable.z(Ids(2));
            obj.display.dendriteMedian(i,4:end) = [X1, Y1, Z1 ];
        obj.display.dendText(i,2) = text(Y1+10, X1+10,num2str(i), ...
            'fontsize',14,'color',[1 1 1],'visible','off');
            V2 = [V2 ; [Y1, X1, Z1]];
        end
        obj.cells.branches(i) = {round(V2./2)};
        continue;
    end
    %% for display
    X1 = obj.cells.dendriteTable.y(Ids);%*2;
    Y1 = obj.cells.dendriteTable.x(Ids);%*2;
    Z1 = obj.cells.dendriteTable.z(Ids);%*2;
    cumDist=sum(sqrt(diff(X1).^2+diff(Y1).^2+ diff(Z1).^2));
    InterpPoints = round(cumDist)*3;
    V=interparc(InterpPoints,X1,Y1,Z1,'linear');
    V2 = unique(round(V),'rows');
    V2Index = sub2ind(size(obj.cells.maskDend),V2(:,1),V2(:,2),V2(:,3));
    Loc1 = round(min(V2,[],1));
    Loc2 = round(max(V2,[],1));
    obj.display.dendriteMedian(i,:) = [Loc1, Loc2];
    obj.display.dendText(i,1) = text(Loc1(2)-10, Loc1(1)-10,num2str(i), ...
        'fontsize',14,'color',[1 1 1],'visible','off');
    obj.display.dendText(i,2) = text(Loc2(2)+10, Loc2(1)+10,num2str(i), ...
        'fontsize',14,'color',[1 1 1],'visible','off');
    obj.cells.maskDend(V2Index) = 1;
    %% for dend center by branch
    X1 = obj.cells.dendriteTable.x(Ids);
    Y1 = obj.cells.dendriteTable.y(Ids);
    Z1 = obj.cells.dendriteTable.z(Ids);
    %%
    cumDist=sum(sqrt(diff(X1).^2+diff(Y1).^2+ diff(Z1).^2));
    InterpPoints = round(cumDist)*3;
    V=interparc(InterpPoints,X1,Y1,Z1,'linear');
    [~, index, ~] = unique(round(V),'rows');
    V2 = round(V(sort(index),:)./2);
    obj.cells.branches(i) = {V2};
end
%%
obj.cells.maskDend = imdilate(obj.cells.maskDend,strel('arbitrary',ones(4,4,4)));
set(obj.handles.ustatH,'String','Ready')
