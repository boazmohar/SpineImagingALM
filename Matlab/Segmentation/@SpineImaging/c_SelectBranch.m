function c_SelectBranch(obj)
Branches = floor(obj.OptionsStruct.SelectedBranches);
obj.All.x = [];
obj.All.y = [];
obj.All.z = [];
obj.All.Id = [];
obj.AllDist=0;
axes(obj.AxesHandle)
cla
hold on;
for i = 1:length(Branches)
    Color = obj.BranchColors(Branches(i),:);
    %% getting which sub section of the branch is wanted
    
    txt = num2str(obj.OptionsStruct.SelectedBranches(i));
    dot = strfind(txt, '.');
    if ~isempty(dot) && isempty(obj.OptionsStruct.MaxBranchDist)
        throw(MException('SpineImaging:c_SelectBranch', ...
            'Do not use dot notation if MaxBranchDist is off'))
    end
    if ~isempty(dot)
        txt = txt(dot+1:end);
        SubBranch = zeros(length(txt),1);
        for s = 1:length(txt)
            SubBranch(s) = str2double(txt(s));
        end
    else
        SubBranch = [];
    end
    %% checking that a subsection selected are in the range possible
    if ~isempty(SubBranch) && max(SubBranch) > length(obj.SubBranch{Branches(i)})+1
        throw(MException('SpineImaging:c_SelectBranch', ...
            sprintf('subselection %d is to large for branch: %d',...
            max(SubBranch),Branches(i))))
    end
    %% plotting all sub brancheing points
    X1 = obj.Centers{Branches(i)}.x;
    Y1 = obj.Centers{Branches(i)}.y;
    Z1 = obj.Centers{Branches(i)}.z;
    if obj.OptionsStruct.PlotReconstraction && ~isempty(obj.OptionsStruct.MaxBranchDist)
        divX = X1(obj.SubBranch{Branches(i)});
        divY = Y1(obj.SubBranch{Branches(i)});
        divZ = Z1(obj.SubBranch{Branches(i)});
        scatter3(divX,divY,-divZ,'s','filled','sizedata',100,...
            'MarkerEdgeColor','b','MarkerFaceColor','b')
        hold on;
    end
    %% adding and padding per sub branch
    if ~isempty(SubBranch)
        CurrentDiv = [0 obj.SubBranch{Branches(i)} length(obj.Centers{Branches(i)}.x)];
        for s = 1:length(SubBranch)
            X2 = X1(CurrentDiv(SubBranch(s))+1:CurrentDiv(SubBranch(s)+1));
            Y2 = Y1(CurrentDiv(SubBranch(s))+1:CurrentDiv(SubBranch(s)+1));
            Z2 = Z1(CurrentDiv(SubBranch(s))+1:CurrentDiv(SubBranch(s)+1));
            if obj.OptionsStruct.PadZ
                [X2, Y2,Z2] = obj.c_PadZ(X2,Y2,Z2);
            elseif isfield(obj.OptionsStruct,'DenseStackMode') && ...
                    obj.OptionsStruct.DenseStackMode
                [X2, Y2,Z2] = obj.c_PadZDense(X2,Y2,Z2);
            end
            obj.All.x = [obj.All.x ; X2];
            obj.All.y = [obj.All.y ; Y2];
            obj.All.z = [obj.All.z ; Z2];
            obj.All.Id = [obj.All.Id ; ones(length(X2),1)*Branches(i)+(SubBranch(s)/10)]; 
            subBranchNum = length(obj.SubBranch{Branches(i)})+1;
            obj.AllDist = obj.AllDist + obj.Dist(Branches(i))/subBranchNum;
            h=scatter3(X2,Y2,-Z2);
            h.CData = Color;
            hold on;
            text(mean(X2)+2,mean(Y2)+2,-mean(Z2),['#' num2str(Branches(i)) '.' ...
                num2str(SubBranch(s))],'FontSize',14, ...
                'Color',Color);
        end
    else
        if isfield(obj.OptionsStruct,'DenseStackMode') && ...
                obj.OptionsStruct.DenseStackMode
            [X1, Y1,Z1] = obj.c_PadZDense(X1,Y1,Z1);
        end
        obj.All.x = [obj.All.x ; X1];
        obj.All.y = [obj.All.y ; Y1];
        obj.All.z = [obj.All.z ; Z1];
        obj.All.Id = [obj.All.Id ; ones(length(X1),1)*Branches(i)]; 
        obj.AllDist = obj.AllDist + obj.Dist(Branches(i));
        h=scatter3(X1,Y1,-Z1);
        h.CData = Color;
        hold on;
        text(mean(X1)+2,mean(Y1)+2,-mean(Z1),['#' num2str(Branches(i))],'FontSize',14, ...
            'Color',Color);
    end
end