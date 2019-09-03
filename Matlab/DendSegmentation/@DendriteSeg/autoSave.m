function autoSave(obj,~,~)
% auto saving using a timer object every 60 sec
if obj.saving.changed > 0
    
    set(obj.handles.ustatH,'String','Auto save')
    pause(0.1)
    if ~exist('c:\dendSave','dir')
        mkdir('c:\dendSave')
    end
    Name = datestr( datetime('now'),'yyyy-mm-dd_HH-MM-SS');
    cells = obj.cells;
    save(['c:\dendSave\' Name '_DendSeg.mat'],'cells','-V7.3')
    disp('auto saved')
    obj.saving.changed = 0;
end

set(obj.handles.ustatH,'String','Ready')
end

