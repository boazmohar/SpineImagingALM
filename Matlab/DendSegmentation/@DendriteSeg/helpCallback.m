function helpCallback(obj,figH,eventData)
%%
Opt.Interpreter = 'tex';
Opt.WindowStyle = 'noraml';
fontName = 'FixedWidth';
fontSize = 14;

msgHandle = msgbox({...
   'Use \bfScrool\rm to move up and down the stack' ...
   'Use \bfLeft click\rm to add a new mask' ...
   'Use \bfShift+Left click\rm to remove a mask',...
   'Use \bfCtrl+Left click\rm to add a border',...
   'Use \bfa\rm to increase contrast', ...
   'Use \bfs\rm to decrease contrast', ...
   'Use \bfz\rm to increase quality', ...
   'Use \bfx\rm to decrease quality', ...
   'Use \bfq\rm to increase contrast ratio', ...
   'Use \bfw\rm to decrease contrast ratio', ...
   'Use \bfe\rm to increase dendrite number', ...
   'Use \bfr\rm to decrease dendrite number', ...
   'Use \bfCtrl+d\rm to enable zoom' ...
   'Use \bfCtrl+f\rm to disable zoom' ...
   'Use \bfCtrl+c\rm to enable pan' ...
   'Use \bfCtrl+v\rm to disable pan' ...
   },'Help', 'none', Opt);


set( msgHandle, 'Visible', 'off' );
% get handles to the UIControls ([OK] PushButton) and Text
kids0 = findobj( msgHandle, 'Type', 'UIControl' );
kids1 = findobj( msgHandle, 'Type', 'Text' );

% change the font and fontsize
extent0 = get( kids1, 'Extent' ); % text extent in old font
set( [kids0, kids1], 'FontName', fontName, 'FontSize', fontSize );
extent1 = get( kids1, 'Extent' ); % text extent in new font

% need to resize the msgbox object to accommodate new FontName
% and FontSize
delta = extent1 - extent0; % change in extent
pos = get( msgHandle, 'Position' ); % msgbox current position
pos = pos + delta*2; % change size of msgbox
set( msgHandle, 'Position', pos ); % set new position

set( msgHandle, 'Visible', 'on' );