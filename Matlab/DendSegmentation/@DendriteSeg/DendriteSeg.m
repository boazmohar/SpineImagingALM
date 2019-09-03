classdef DendriteSeg < handle
    %DendriteSeg GUI for segmentation of spines
    %   give a 3d array for the image as the first argument or you will be
    %   promted to select a tif file. Will auto-save to c:\dendSave every
    %   minute. press Ctrl+h for help with keyboard shortcuts
    
    properties
        handles = struct() % for all graphics info
        cells = struct() % info to be saved about the spines 
        display = struct() % image and othe information
        current = struct() % current spine to add
        saving = struct() % info related to saving and loading
    end
    
    methods
        function obj = DendriteSeg(varargin)
            %% init some variables
            obj.display.pixelSize = [0.2 0.2 0.5];
            obj.display.Im_contrast = 1;
            obj.display.Z = 3;
            obj.display.showDend = 0;
            obj.initFigure();       
            %% set auto save every minute
            obj.saving.changed = 0;
            delete(timerfind('Tag', 'DendriteSeg timer'));
            obj.saving.t = timer('Tag', 'DendriteSeg timer', 'Period', 720,...
                'ExecutionMode', 'fixedDelay', 'StartDelay', 720, ...
                'TimerFcn', @obj.autoSave, 'TasksToExecute', inf);
            start(obj.saving.t);
            
            %% load image file or have it as  the first input
            if nargin > 0
                obj.display.volume = varargin{1};
            else
                obj.loadImage();
            end
            %% rotate image to better fit the screen
%             if size(obj.display.volume,1) > size(obj.display.volume,2)
%                 obj.display.volume = permute(obj.display.volume,[2,1,3]);
%             end
            %% init graphics and drew first image
            obj.setupImage()
            %% init cells
            obj.cells.labelSpine = zeros(size(obj.display.volume));
            obj.cells.nTotal = 0;
            set(obj.handles.ustatH,'String','Ready')
        end
        %% function to delete the timer when the object is destroyd
        function delete(obj)
            delete(obj.saving.t);
        end
        %% init functions
        loadImage(obj)
        initFigure(obj)
        setupImage(obj)
        %% update functions
        updateImage(obj)
        imgOut = normImage(obj,imgIn)
        imShade(obj)
        %% callbacks
        buttonCallback(obj,figH,eventData)
        keyCallback(obj,figH,eventData)
        scrollCallback(obj,figH,eventData)
        helpCallback(obj,figH,eventData)
        loadCallback(obj,figH,eventData)
        masksCallback(obj,figH,eventData)
        saveCallback(obj,figH,eventData)
        autoSave(obj,timerObj,event)
        %% dendrite mask functions
        swcCallback(obj,figH,eventData)
        dendReadSWC(obj)
        dendAddBranches(obj)
        dendMaskUpdate(obj)
        %% processing
        Error = checkLimits(obj,X,Y)
        Error = getMeta(obj)
        addCell3d(obj)
        [bwNew, bwCell, cellMean] = subAddCell(obj, bwCurr, fovImg, cellMean)
        addBorder(obj)
        deleteMask(obj,X,Y,Z)
    end
end

