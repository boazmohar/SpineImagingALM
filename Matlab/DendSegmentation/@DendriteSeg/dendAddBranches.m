function dendAddBranches(obj)
set(obj.handles.ustatH,'String','AddBranches')
Base = find(obj.cells.dendriteTable.pNum == -1);
Branchs = zeros(height(obj.cells.dendriteTable),1);
obj.cells.dendriteTable.BranchID = Branchs;
BranchID = 0;
%% go over the tree
while ~isempty(Base)
    BranchID = BranchID+1;
    obj.cells.dendriteTable.BranchID(Base) = BranchID;
    Parents = find(obj.cells.dendriteTable.pNum == Base(1));
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
    obj.cells.dendriteTable.BranchID(Current) = BranchID;
    Flag = 1;
    while Flag
        Parents = find(obj.cells.dendriteTable.pNum == Current);
        if isempty(Parents)
            Flag = 0;
        elseif length(Parents) > 1
            Base = [Base; Parents(2:end)];
            Current = Parents(1);
            BranchID = BranchID+1;
            obj.cells.dendriteTable.BranchID(Current) = BranchID;
        else
            Current = Parents(1);
            obj.cells.dendriteTable.BranchID(Current) = BranchID;
        end
    end
end
obj.cells.BranchNum = BranchID;