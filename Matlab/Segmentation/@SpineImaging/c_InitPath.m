function c_InitPath(obj)
%  initPath(obj) updates the objects centers
obj.AllDist =0;
obj.All.x = [];
obj.All.y = [];
obj.All.z = [];
obj.All.Id = [];
axes(obj.AxesHandle)
cla;
hold on;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
for i = 1:obj.BranchNum
    Color = obj.BranchColors(i,:);
    %% interpolation and smoothing
    Ids = find(obj.Table.BranchID == i);
    X1 = obj.Table.x(Ids);
    Y1 = obj.Table.y(Ids);
    Z1 = obj.Table.z(Ids);
    
    % add the privous one if possible
    PrevId = obj.Table.pNum(Ids(1));
    if PrevId ~= -1
        X1 = [obj.Table.x(PrevId); X1];
        Y1 = [obj.Table.y(PrevId); Y1];
        Z1 = [obj.Table.z(PrevId); Z1];
    end
    if obj.OptionsStruct.PlotReconstraction
        h = scatter3(X1,Y1,-Z1);
        h.CData = Color;
        text(mean(X1)+2,mean(Y1)+2,-mean(Z1),['#' num2str(i)],'FontSize',14, ...
            'Color',Color);
    end
    cumDist=sum(sqrt(diff(X1).^2+diff(Y1).^2+ diff(Z1).^2));
    InterpPoints = round(cumDist./obj.OptionsStruct.InterpPointsDist);
    if length(X1) > 1 % if branch is bigger then one point
        V=interparc(InterpPoints,X1,Y1,Z1,'linear');
        LowpassWindow = round(cumDist./obj.OptionsStruct.LowpassWindowUm);
        LowpassWindow = max([min([LowpassWindow+1,length(V(:,1))-2]) 3]);
        LowpassWindow = 2*floor(LowpassWindow/2)+1;
        %LowpassWindow= min([LowpassWindow, length(V(:,1))]);
        X=sgolayfilt(V(:,1),1,LowpassWindow);
        Y=sgolayfilt(V(:,2),1,LowpassWindow);
        Z=sgolayfilt(V(:,3),1,LowpassWindow);
    else % else take the point as is
        X=X1;  Y=Y1; Z = Z1;
    end
    if obj.OptionsStruct.PlotLine
        p = plot3(X,Y,-Z);
        if ~isempty(Color)
            p.Color = Color;
        else
            Color = p.Color;
        end
    end
    %% distence function in 3d assuming square X=Y FOV
    if length(X) > 3
        Zprime = Z.*obj.OptionsStruct.zAspectRatio; % weighing Z to ensure spacing with vertical dendrites
        Xprime = X./(obj.OptionsStruct.xSize ./obj.OptionsStruct.ySize);
        dist=sqrt(diff(Xprime).^2+diff(Y).^2+ diff(Zprime).^2);
        %% find center of volumes
        cDist = cumsum(dist);
        Vols = rem(cDist+obj.OptionsStruct.StepSize/2,obj.OptionsStruct.StepSize);
        [~,VolLocs] = findpeaks(Vols);
        CenterX = X(VolLocs);
        CenterY = Y(VolLocs);
        CenterZ = Z(VolLocs);
        %% subdividing into sub branches
        if ~isempty(obj.OptionsStruct.MaxBranchDist) && cumDist > obj.OptionsStruct.MaxBranchDist
            subBranchNum = ceil(cumDist/ obj.OptionsStruct.MaxBranchDist);
            divPoint = round(length(CenterX)/subBranchNum);
            for s = 2:subBranchNum-1
                divPoint(s) = divPoint(s-1)+ round(length(CenterX)/subBranchNum);
            end
            obj.SubBranch{i} =divPoint;
            if obj.OptionsStruct.PlotReconstraction
                divX = CenterX(obj.SubBranch{i});
                divY = CenterY(obj.SubBranch{i});
                divZ = CenterZ(obj.SubBranch{i});
                scatter3(divX,divY,-divZ,'s','filled','sizedata',100,...
                    'MarkerEdgeColor','b','MarkerFaceColor','b')
            end
        else
             obj.SubBranch{i} = [];
        end
    else
        cumDist = 0;
        CenterX = X;
        CenterY = Y;
        CenterZ = Z;
    end
    %% if pad z
    if obj.OptionsStruct.PadZ && isempty(obj.OptionsStruct.SelectedBranches)
        [CenterX, CenterY,CenterZ] = obj.c_PadZ(CenterX,CenterY,CenterZ);
    end
    %%
    if obj.OptionsStruct.PlotCenters
        h2=scatter3(CenterX,CenterY,-CenterZ,24,'filled');
        if ~isempty(Color)
            h2.CData = Color;
        end
    end
    %% add to global
    obj.All.x = [obj.All.x; CenterX];
    obj.All.y = [obj.All.y; CenterY];
    obj.All.z = [obj.All.z; CenterZ];
    obj.All.Id = [obj.All.Id ; ones(length(CenterX),1)*i]; 
    obj.Centers{i}.x = CenterX;
    obj.Centers{i}.y = CenterY;
    obj.Centers{i}.z = CenterZ;
    obj.Dist(i) = cumDist;
    obj.AllDist = obj.AllDist +cumDist;
end