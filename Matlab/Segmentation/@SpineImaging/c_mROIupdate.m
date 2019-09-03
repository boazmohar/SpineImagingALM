function c_mROIupdate( obj,Options )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% update options strct
try
    axes(obj.AxesHandle)
    obj.TextHandle.Visible = 'On';
    obj.TextHandle.String = 'Working';
    obj.OptionsStruct = Options;
    obj.linesPerField = obj.OptionsStruct.ImagingLines + obj.OptionsStruct.FlyLines;
    %if I change the file nothing happens!!!
    %% Get centers by branch
%     if ~obj.OptionsStruct.BrushingMode
    obj.c_InitPath();
%     end
    %% select braches?
    if ~isempty(obj.OptionsStruct.SelectedBranches)
        obj.c_SelectBranch();
    elseif isfield(obj.OptionsStruct,'DenseStackMode') && ...
            obj.OptionsStruct.DenseStackMode
        [obj.All.x, obj.All.y,obj.All.z] = ...
        obj.c_PadZDense(obj.All.x,obj.All.y,obj.All.z);
    end
    %% Brushing mode
%     if obj.OptionsStruct.BrushingMode
%         axes(obj.AxesHandle)
%         cla
%         hold on;
%         S = scatter3(obj.All.x,obj.All.y,-obj.All.z,'filled');
%         obj.FigureHandle.Units = 'Pixels';
%         obj.AxesHandle.Units =  'Pixels';
%         brush on;
%         obj.BrushingMask = S.BrushData;
%     else
%         brush off;
%         rotate3d on;
%     end
    %% Z Re-profile
    if obj.OptionsStruct.ZOptimize
        obj.c_ZOptimize();
    end
    %% Update display
    ScanFields = length(obj.All.x);
    framerate = 1/(obj.linesPerField/obj.OptionsStruct.lineRate*ScanFields);
    TitleText = ['Scanfields: ' num2str(ScanFields) ', Length: ' num2str(obj.AllDist,3) ...
        'um, Vol rate: ' num2str(framerate,3) 'Hz'];
    Zspan = max(obj.All.z) - min(obj.All.z);
    axes(obj.AxesHandle)
    if Zspan > obj.OptionsStruct.ZmaxTravel
        TitleText = [TitleText ', Zspan > max: ' num2str(Zspan,4) 'um'];
    else
        TitleText = [TitleText ', Zspan: ' num2str(Zspan,4) 'um'];
    end
    if ~isempty(obj.maxErrorMicron)
        TitleText = [TitleText ', maxErr: ' num2str(obj.maxErrorMicron,3)];
    end
    title(TitleText);
    obj.TextHandle.Visible = 'Off';
catch E
    obj.TextHandle.String = ['Error: ' E.message];
    rethrow(E);
end