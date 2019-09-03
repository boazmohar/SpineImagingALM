function InitFigure(obj)
%% target figure
obj.TargetHandle = figure();
obj.TargetHandle.Units = 'Normalized';
obj.TargetHandle.Position = [0.05 0.05 0.9 0.85];
set(obj.TargetHandle,'name','TARGET VIEW','numbertitle','off');
obj.TargetHandle.Units = 'pixels';
colormap(obj.TargetHandle,'gray');
cameratoolbar();
obj.TargetAxesHandle = axes('Parent',obj.TargetHandle,'Color',[0 0 0 ]);
view(obj.TargetAxesHandle,37,45);
U1 = uicontrol('Style','slider','Callback',@changeAlpha, ...
        'Position', [20 940 150 50], 'Min',-Inf,'Max',Inf,'Value',0,...
        'String','Low');
U2 = uicontrol('Style','slider','Callback',@changeAlpha, ...
        'Position', [20 840 150 50], 'Min',-Inf,'Max',Inf,'Value',1,...
        'String','High');
obj.TargetSliders(1) =U1;
obj.TargetSliders(2) =U2;
    
%% main figure
obj.FigureHandle = figure();
obj.FigureHandle.Units = 'Normalized';
obj.FigureHandle.Position = [0.05 0.05 0.9 0.85];
obj.FigureHandle.Color = [0.9 0.9 0.9];
obj.PanelHandle = uipanel('parent',obj.FigureHandle','position',[0 0 0.3 1]);

set(obj.FigureHandle,'name','SPINE IMAGEING','numbertitle','off');
uicontrol('parent',  obj.PanelHandle,'position',[30 7 100 22], ...
    'Style','pushbutton','Callback',@obj.c_AddFields,'string','Add ROIs');
uicontrol('parent', obj.PanelHandle,'position',[150 7 100 22], ...
    'Style','pushbutton','Callback',@obj.c_SaveObjCallback,'string','Save');
uicontrol('parent', obj.PanelHandle,'position',[270 7 100 22], ...
    'Style','pushbutton','Callback',@obj.c_displayTarget,'string','Target');
obj.AxesHandle = axes('parent',obj.FigureHandle,'position',[0.35 0.05 0.61 0.9]);
obj.AxesHandle.Color = [0.9 0.9 0.9];
obj.Axes2Handle =axes('parent', obj.PanelHandle,'position',[0.09 0.07 0.8 0.27]);
obj.TextHandle = uicontrol('style','text','units','normalized', ...
    'position',[0.3,0.89 0.1 0.1],'fontsize',14','foregroundcolor',...
    [1 0 0],'String','Working','parent', obj.FigureHandle);
obj.FigureHandle.Units = 'pixels';
obj.BranchColors = distinguishable_colors(200);
set(obj.AxesHandle,'xlimmode','manual','ylimmode','manual')
set(obj.AxesHandle,'xlim',[0 530],'ylim',[0 530])
axes(obj.AxesHandle);
view(3);
end
function changeAlpha(hObject,~)
   	obj = evalin('base','Sp');
    newval =hObject.Value;
    Alim = obj.TargetAxesHandle.ALim;
    if strcmp(hObject.String,'Low')
        Alim(1) = newval;
    else
        Alim(2) = newval;
    end
    obj.TargetAxesHandle.ALim = Alim;
end
