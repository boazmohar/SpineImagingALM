classdef SpineImaging < handle
    %SpineImaging is an interface to control mROI mode in ScanImage from a SWC
    %reconstraction file
    properties
        %% Options
        OptionsStruct       % sets options for all the classes functions
        OptionsSpec         % sets how these options will be displyed
        OptionsObj          % an object that creates the liseners for the changes in values
        % and calles c_mROIupdate
        %% Data
        Table               % stores x y z and perent data from the SWC file
        Centers             % a cell array of sturcters (x,y,z) by branch number
        % of interpulated centers of the dendrite
        SubBranch           % a cell array of indexes in Centers to sub divide each
        % branch that is more than MaxBranchDist
        BranchColors        % set of colors for each branch
        All                 % a struct (x,y,z) of all the current centers to scan
        Roi                 % a struct (x,y,size) of the normalized ROIs positions
        Dist                % Distence by b ranc
        AllDist             % the current length of dnedrite to scan
        BranchNum           % number of branches detected in the SWC file
        BrushingMask        % (not used) for data brushing mode
        Powers              % a laser power value for each field
        Stack               % the anatomy stack
        %% Graphics
        FigureHandle
        TargetHandle
        TargetAxesHandle
        TargetSliders
        PanelHandle
        AxesHandle
        Axes2Handle
        TextHandle
        %% Z optimization
        Zwave
        ZError
        TrueFieldsMask
        TrueFieldsCenterSamp%
        CmdWave             % output Z to SI
        maxErrorMicron
        linesPerField       % computed from options lines to scan and overhead
        Offset              % offset in microns in Z due to field centering
    end
    %% constractor function
    methods
        function obj = SpineImaging(varargin)
            obj.c_SetDefaultOptions();        % create the defulat option
            obj.InitFigure();           % initialize grahics
            switch nargin
                case 0 % no input arguments request a SWC file
                    [Name,Path] = uigetfile('*.SWC');
                    obj.OptionsStruct.filename = [Path Name];
                case 1 % first input is the filename to open
                    obj.OptionsStruct.filename = varargin{1};
                case 2 % from loading
                    s = varargin{2};
                    fields = fieldnames(s);
                    for i=1:numel(fields)
                        if isempty(strfind(fields{i},'Handle')) && ...
                                ~strcmpi(fields{i}, 'OptionsObj')&& ...
                                ~strcmpi(fields{i}, 'TargetSliders');
                            obj.(fields{i}) = s.(fields{i});
                        end
                    end
                    if isfield(s.OptionsStruct,'MaxBranchDist')
                        obj.OptionsStruct.MaxBranchDist = s.OptionsStruct.MaxBranchDist;
                    else
                        obj.OptionsStruct.MaxBranchDist = [];
                    end
                otherwise
                    throw(MException('SpineImaging:Constractor',...
                        'Number of variables (%d) wrong',nargin));
            end
            if nargin ~= 2
                
                obj.c_ReadSWC();            % read SWC file into Table
                obj.c_DetectBranches();     % detect branches in the file
            end
            obj.OptionsObj = fn_control(obj.OptionsStruct,@obj.c_mROIupdate,...
                obj.OptionsSpec,obj.PanelHandle);   % create the options obj
            c_mROIupdate(obj,obj.OptionsStruct);    % do the first update manually
            Zlim = -[max(obj.Table.z)*1.1 min(obj.Table.z)*0.8 ];
            set(obj.AxesHandle,'zlimmode','manual','zlim',Zlim)
            obj.TextHandle.Visible = 'Off';
        end
        %% save and load functions
        function  s = saveobj(obj)
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
            s = struct(obj);
            obj.AxesHandle          = A;
            obj.PanelHandle         = P;
            obj.FigureHandle        = F;
            obj.OptionsObj          = O;
            obj.TargetHandle        = TF;
            obj.TargetAxesHandle    = TA;
            obj.TargetSliders       = TS;
        end
        c_ReadSWC(obj)
        c_DetectBranches(obj)
    end
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                obj = SpineImaging('load',s);
            else
                obj = s;
            end
        end
    end
    %% Other methods
    methods (Access = 'private')
        [X, Y, Z] = c_PadZDense(obj,X,Y,Z)
        c_SelectBranch(obj)
        c_SetDefaultOptions(obj)
        c_mROIupdate(obj,Options)
        traj=c_ZOptAna(obj,Amax,Vmax,SampRate,V0,V1,dist,Verbose)
        c_ZOptimize(obj)
        c_ConvertROIs(obj);
        c_PiezoOptIterative(obj);
        c_AddROIs(obj);
        c_SaveObjCallback(obj,hObject, eventdata, handles);                         % do the first update manually
        c_AddFields(obj,hObject, eventdata, handles)
        c_CopyProps(obj,s);
        [CenterX, CenterY,CenterZ] = c_PadZ(obj,CenterX,CenterY,CenterZ);
        c_SetPower(obj)
        c_displayTarget(obj,hObject, eventdata, handles)
    end
end

