%% load Sp.mat
obj = Sp_s;
fig = figure();
clf()
branches =obj.OptionsStruct.SelectedBranches;
UM_Zoom = obj.OptionsStruct.UM_1X./obj.OptionsStruct.Zoom;
Pixel_Size = UM_Zoom./obj.OptionsStruct.RefPixels;

for i = 1:obj.BranchNum
    Ids = find(obj.Table.BranchID == i);
    X1 = obj.Table.x(Ids)./Pixel_Size;
    Y1 = obj.Table.y(Ids)./Pixel_Size;
    Z1 = obj.Table.z(Ids)./obj.OptionsStruct.Zstep;
    
    % add the privous one if possible
    PrevId = obj.Table.pNum(Ids(1));
    if PrevId ~= -1
        X1 = [obj.Table.x(PrevId)./Pixel_Size; X1];
        Y1 = [obj.Table.y(PrevId)./Pixel_Size; Y1];
        Z1 = [obj.Table.z(PrevId)./obj.OptionsStruct.Zstep; Z1];
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
    p = plot3(X,Y,-Z, 'color','white');
    hold('on')
    if sum(i==branches)
        p.Color = 'red';
        p.LineWidth=3;
        
    else
        p.Color = [0.7 0.7 0.7];
    end
    if i == 1
        scatter3(X, Y, -Z, 50, 'b', 'filled')
    end
end
%
xlabel('Y (um)')
ylabel('X (um)')
zlabel('Z (um)')

%
view(2)
