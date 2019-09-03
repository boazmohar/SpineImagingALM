function swcCallback(obj,figH,eventData)
if strcmp(get(eventData.Source,'Checked'),'off')
    set(eventData.Source,'Checked','on')
    obj.display.showDend = 1;
else
    set(eventData.Source,'Checked','off')
    obj.display.showDend = 0;
end
obj.updateImage();