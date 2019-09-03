function c_AddFields(obj,hObject, eventdata, handles)
try
    % make sure we don't have two exect same Z positions
    axes(obj.AxesHandle)
    obj.TextHandle.Visible = 'On';
    obj.TextHandle.String = 'Working';
    drawnow
    pause(0.1);
    %% Convert Centers to SI mROI = [x y width height]. Units are scan fraction (0-1) of the reference field of view
    obj.c_ConvertROIs();
    %% Optimizting Z piezo response to hit target positions
    
    obj.c_PiezoOptIterative();
    %% setting power for all scanfields
%     obj.c_SetPower()
    %% Getting SI ready and adding ROIs (add Z waveform setup)
    obj.c_AddROIs();
   
    %% set title
    axes(obj.AxesHandle)
    Title = get(gca,'title');
    index = strfind(Title.String,', maxErr: ');
    if ~isempty(index)
        Title.String = Title.String(1:index-1);
    end
    Title.String= [Title.String ', maxErr: ' num2str(obj.maxErrorMicron,2) ...
        ', Offset: ' num2str(obj.Offset,4)];
    set(gca,'Title',Title)
    %% automatic saving
    try
        cd('D:\Data\SpSave')
        Name = datestr( datetime('now'),'yyyy-mm-dd_HH-MM-SS');
        [~,file,~] = fileparts(obj.OptionsStruct.filename);
        Name = [Name '_' file];
        % save handels
        A   = obj.AxesHandle;
        P   = obj.PanelHandle;
        F   = obj.FigureHandle;
        O   = obj.OptionsObj;
        TF  = obj.TargetHandle;
        TA  = obj.TargetAxesHandle;
        TS  = obj.TargetSliders;
        % delete stack
        if isprop(obj,'Stack')
            ST = obj.Stack;
            obj.Stack = [];
        end
        % convert to struct
        obj.AxesHandle          = [];
        obj.PanelHandle         = [];
        obj.FigureHandle        = [];
        obj.OptionsObj          = [];
        obj.TargetHandle        = [];
        obj.TargetAxesHandle    = [];
        obj.TargetSliders       = [];
        Sp_s = struct(obj);
        save(Name,'Sp_s');
        obj.AxesHandle          = A;
        obj.PanelHandle         = P;
        obj.FigureHandle        = F;
        obj.OptionsObj          = O;
        obj.TargetHandle        = TF;
        obj.TargetAxesHandle    = TA;
        obj.TargetSliders       = TS;
        if isprop(obj,'Stack')
            obj.Stack           = ST;
        end
    catch e
        disp(e.message);
    end
     obj.TextHandle.Visible = 'Off';
catch E
    obj.TextHandle.String = 'Error';
    rethrow(E);
end