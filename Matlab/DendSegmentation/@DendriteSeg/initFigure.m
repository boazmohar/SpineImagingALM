function initFigure(obj)
%% open figure 
figH = figure('units','normalized','position',[0.15 0.15 0.7 0.7],'menubar','none');
axH = axes('Position', [0.05 0.15 0.90 0.84], 'parent',figH);
set(axH,'units','pixels')
set(figH,'units','pixels')
figColor = get(figH, 'Color');
set(axH, 'Units', 'normalized');
%% Contrast threshold
uicontrol('Style', 'text', ...
    'String', {'Contrast threshold:', ' fract of mean'}, ...
    'Units', 'pixels', 'Position', [5 45, 80, 40], ...
    'BackgroundColor', figColor, 'HorizontalAlignment', 'left');
utH = uicontrol('Style', 'edit', 'String', '0.95',  'Units', 'pixels', ...
    'Position', [5 5, 60, 30]);
%% Dilation disk
uicontrol('Style', 'text', 'String', {'Dilation disk:', ' size in pix'}, ...
    'Position', [85 45, 80, 30], 'BackgroundColor', figColor, ...
    'HorizontalAlignment', 'left');
udH = uicontrol('Style', 'edit',  'String', '1', ...    % default diskR
    'Units', 'pixels', 'Position', [85 5, 60, 30]);
%% Dendrite Number
uicontrol('Style', 'text', 'String', 'Dend #', 'Position', [170, 35, 80, 30], ...
    'BackgroundColor', figColor, 'HorizontalAlignment', 'left');
unH = uicontrol('Style', 'edit', 'String', '1', ...    % default dendrite number
    'Units', 'pixels', 'Position', [170, 5, 60, 30]);
%% Quality
uicontrol('Style', 'text', 'String', 'Quality', 'Position', [255 35, 80, 30], ...
    'BackgroundColor', figColor, 'HorizontalAlignment', 'left');
uqH = uicontrol('Style', 'popup', 'String', {'High','Mid','Low','Other'}, ...
    'Units', 'pixels',  'Position', [255 5, 60, 30]);
%% Notes
uicontrol('Style', 'text', 'String', 'Notes', 'Position', [330, 35, 80, 30], ...
    'BackgroundColor', figColor, 'HorizontalAlignment', 'left');
unoteH = uicontrol('Style', 'edit', 'String', '', ...    % default dendrite number
    'Units', 'pixels', 'Position', [330 5, 400, 40]);
%% Status
uicontrol('Style', 'text', 'String', 'Status', 'Position', [750, 35, 80, 30], ...
    'BackgroundColor', figColor, 'HorizontalAlignment', 'left','fontsize',14);
ustatH = uicontrol('Style', 'text', 'String', 'Init', ...    % default dendrite number
    'Units', 'pixels', 'Position', [750 5, 80, 40],'ForegroundColor','r','fontsize',14);
%% Z position
uicontrol('Style', 'text', 'String', 'Z', 'Position', [1060, 25, 50, 30], ...
    'BackgroundColor', figColor, 'HorizontalAlignment', 'center','fontsize',14);
uzH = uicontrol('Style', 'text', 'String', num2str(obj.display.Z), ...    % current Z
    'Units', 'pixels', 'Position', [1110,25, 50, 30], 'HorizontalAlignment', 'center',...
    'fontsize',14);
%% save handels
obj.handles.figH = figH;
obj.handles.axH = axH;
obj.handles.figColor = figColor;
obj.handles.utH = utH;
obj.handles.udH = udH;
obj.handles.unH = unH;
obj.handles.uqH = uqH;
obj.handles.unoteH = unoteH;
obj.handles.ustatH = ustatH;
obj.handles.uzH = uzH;
%% file menu
mymenu2 = uimenu('Parent',figH,'Label','File');
uimenu('Parent',mymenu2,'Label','Load mask', 'callback', @obj.loadCallback);
uimenu('Parent',mymenu2,'Label','Save mask', 'callback', @obj.saveCallback);
uimenu('Parent',mymenu2,'Label','Load dendrite', 'callback', @obj.swcCallback);
uimenu('Parent',mymenu2,'Label','Show dendrite', 'callback', @obj.dendCallback);
uimenu('Parent',mymenu2,'Label','Prepare Masks', 'callback', @obj.masksCallback);

%% hotkeys menu
mymenu = uimenu('Parent',figH,'Label','Hot Keys');
uimenu('Parent',mymenu,'Label','Zoom on','Accelerator','d','Callback',@(src,evt)zoom(figH,'on'));
uimenu('Parent',mymenu,'Label','Zoom off','Accelerator','f','Callback',@(src,evt)zoom(figH,'off'));
uimenu('Parent',mymenu,'Label','Pan','Accelerator','c','Callback',@(src,evt)pan(figH,'on'));
uimenu('Parent',mymenu,'Label','Pan off','Accelerator','v','Callback',@(src,evt)pan(figH,'off'));
%% help menu
mymenu3 = uimenu('Parent',figH,'Label','Help');
uimenu('Parent',mymenu3,'Label','Show help','Accelerator','h','Callback',@obj.helpCallback);
%% callbacks
set(figH, 'WindowButtonDownFcn', @obj.buttonCallback, ...
             'KeyPressFcn', @obj.keyCallback, ...
             'WindowScrollWheelFcn', @obj.scrollCallback);
%% params fn_control

end