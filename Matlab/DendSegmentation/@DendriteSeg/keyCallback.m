function keyCallback(obj,~,eventData)
% callback function for key pressing
switch lower(eventData.Key)
    case 'a' % increase brightset
        obj.display.brightset= obj.display.brightset.*1.1;
        obj.updateImage()
    case 's' % decrease brightset
        obj.display.brightset= obj.display.brightset.*0.9;
        obj.updateImage()
    case 'z' % move down in Z
        % remap to change the q
        value = get(obj.handles.uqH,'value');
        if value > 1
            value = value - 1;
        end
        set(obj.handles.uqH,'value', value);
    case 'x' % move up in Z
        value = get(obj.handles.uqH,'value');
        if value < 4
            value = value + 1;
        end
        set(obj.handles.uqH,'value', value);
    case 'q' %decrease threshold
        cThresh = str2double(get(obj.handles.utH, 'String'));
        cThresh = cThresh-0.05;
        if ~isnan(cThresh)
            set(obj.handles.utH, 'String',num2str(cThresh))
        end
    case 'w' % increase threshold
        cThresh = str2double(get(obj.handles.utH, 'String'));
        cThresh = cThresh+0.05;
        if ~isnan(cThresh)
            set(obj.handles.utH, 'String',num2str(cThresh))
        end
    case 'e' %decrease threshold
        dend = str2double(get(obj.handles.unH, 'String'));
        dend = dend-1;
        dend = round(dend);
        if isfinite(dend) && dend > 0
            set(obj.handles.unH, 'String',num2str(dend))
        end
    case 'r' % increase threshold
        dend = str2double(get(obj.handles.unH, 'String'));
        dend = dend+1;
        dend = round(dend);
        if isfinite(dend) && dend > 0
            set(obj.handles.unH, 'String',num2str(dend))
        end
end
end