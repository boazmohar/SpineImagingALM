function c_DetectBranches( obj )
% DetectBranches( obj ) adds a branch id number to the obhects table
%% find first base segment
Base = find(obj.Table.pNum == -1);
% if isempty(Base) || length(Base) > 1
%     throw(MException('SpineImaging:c_DetectBranches','Error with the tree: check disconnected branches'));
% end
%% init branch id coulmn of the table
Branchs = zeros(height(obj.Table),1);
obj.Table.BranchID = Branchs;
BranchID = 0;
%% go over the tree
while ~isempty(Base)
    BranchID = BranchID+1;
    obj.Table.BranchID(Base) = BranchID;
    Parents = find(obj.Table.pNum == Base(1));
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
    obj.Table.BranchID(Current) = BranchID;
    Flag = 1;
    while Flag
        Parents = find(obj.Table.pNum == Current);
        if isempty(Parents)
            Flag = 0;
        elseif length(Parents) > 1
            Base = [Base; Parents(2:end)];
            Current = Parents(1);
            BranchID = BranchID+1;
            obj.Table.BranchID(Current) = BranchID;
        else
            Current = Parents(1);
            obj.Table.BranchID(Current) = BranchID;
        end
    end
end
obj.BranchNum = BranchID;