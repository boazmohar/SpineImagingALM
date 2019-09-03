function c_displayTarget(obj,hObject, eventdata, handles)
%displayTarget display a target from the anatomy stack

%% update options strct
xScale=1;
yScale=1;
try
    
    axes( obj.TargetAxesHandle);
    cla;
    obj.TextHandle.Visible = 'On';
    obj.TextHandle.String = 'Working';
    if ~obj.OptionsStruct.ZOptimize
        error('ZOptimize is false')
    end
    if isempty(obj.Stack)
        if isfield(obj.OptionsStruct,'Stack') && ~isempty(obj.OptionsStruct.Stack)
            [pathname, filename] = fileparts(obj.OptionsStruct.Stack);
        else
            [filename,pathname] = uigetfile('*.tif');
            obj.OptionsStruct.Stack = [pathname filename];
            
        end
        stack=readtiff(pathname,[],filename(1:end-4));
        obj.Stack = stack;
    else
        stack = obj.Stack;
    end
    
    xSize=obj.OptionsStruct.xSize*xScale; %xGalvo brings to center of res
    ySize=obj.OptionsStruct.ySize;%yGalvo brings to start of y, wait to scale
    
    %get scanfield centers
    All.x = obj.All.x;
    All.y = obj.All.y;
    All.z = obj.All.z;
    fieldCenters=table2array(struct2table(All));
    fieldCenters=fieldCenters(obj.TrueFieldsMask==1,:);
    nFields=size(fieldCenters,1);
    
    xAnatPixSize=obj.OptionsStruct.UM_1X./obj.OptionsStruct.RefPixels;
    yAnatPixSize=obj.OptionsStruct.UM_1X./obj.OptionsStruct.RefPixels;
    zAnatPixSize=obj.OptionsStruct.Zstep;
    
    fieldCorners=NaN(nFields,3);
    fieldCorners(:,1)=round((fieldCenters(:,1)-xSize./2)./xAnatPixSize);
    fieldCorners(:,2)=round((fieldCenters(:,2)-ySize./2)./yAnatPixSize);
    fieldCorners(:,3)=round(fieldCenters(:,3)./zAnatPixSize);
    
    ySize=ySize*yScale;
    
    xPix=round(xSize./xAnatPixSize);
    yPix=round(ySize./yAnatPixSize);
    
    eqFields=NaN(yPix,xPix,nFields,'single');
    for i=1:nFields
        xOffsets=[1-min([1,fieldCorners(i,1)]) obj.OptionsStruct.RefPixels-max([obj.OptionsStruct.RefPixels,fieldCorners(i,1)+xPix-1])];
        yOffsets=[1-min([1,fieldCorners(i,2)]) obj.OptionsStruct.RefPixels-max([obj.OptionsStruct.RefPixels,fieldCorners(i,2)+yPix-1])];
        eqFields((1+yOffsets(1)):(end+yOffsets(2)),(1+xOffsets(1)):(end+xOffsets(2)),i)=...
            stack(...
            (fieldCorners(i,2)+yOffsets(1)):(fieldCorners(i,2)+yPix-1+yOffsets(2)),...
            (fieldCorners(i,1)+xOffsets(1)):(fieldCorners(i,1)+xPix-1+xOffsets(2)),...
            fieldCorners(i,3)+1);
    end
    
    %scaling
    eqFieldsResize=NaN([obj.OptionsStruct.ImagingLines,obj.OptionsStruct.ImagingPixels,nFields]);
    for i=1:nFields
        eqFieldsResize(:,:,i)=imresize(eqFields(:,:,i),[obj.OptionsStruct.ImagingLines obj.OptionsStruct.ImagingPixels]);
    end
    
    %%
    Mean =eqFieldsResize;
    maxMean = max(Mean(:));
    minMean = min(Mean(:));
    [xSize, ySize, zSize] = size(Mean);
    %
    hAx = obj.TargetAxesHandle;
    %% set up surface for texture mapping
    for i = 1:zSize
        data = squeeze(Mean(:,:,i));
        X = obj.All.x(i);
        Y = obj.All.y(i);
        Z = obj.All.z(i);
        [x,y,z] = meshgrid([X X+xSize],[Y Y+ySize],-Z);
        hSurf = surface('Parent',hAx,'XData',x,'YData',y,'ZData',z,...
            'FaceColor','texturemap','CDataMapping','scaled',...
            'FaceLighting','none', 'EdgeColor','none');
        set(hSurf,'CData',data,'AlphaData',data,'AlphaDataMapping','scaled',...
            'FaceAlpha','texturemap'); % set surface image data
    end
    set(obj.TargetSliders(1),'Min',minMean,'Max',maxMean,'value',minMean);
    set(obj.TargetSliders(2),'Min',minMean,'Max',maxMean,'Value',maxMean);
    obj.TextHandle.Visible = 'Off';
catch E
    obj.TextHandle.String = ['Error: ' E.message];
    obj.OptionsStruct.Stack  = [];
    rethrow(E);
end
end
