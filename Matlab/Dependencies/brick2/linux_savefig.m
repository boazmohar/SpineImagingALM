function linux_savefig(hf,fname)
% function linux_savefig(hf[,fname])
%---
% save figure with handle hf in folder fn_cd('capture'), under the name
% figuretitle_date_time.png
%
% works only under linux, with ImageMagick installed

% Thomas Deneux
% Copyright 2009-2012

if nargin==0, help linux_savefig, return, end



if strfind(computer,'MAC')
    
    % Mac: use 'screencapture' command to capture the full screen
    figure(hf), pause(.1)
    system(['screencapture ''' fname '''']);
    
    % Cut the image according to figure position
    a = imread(fname);
    pos = get(hf,'position'); pos(2) = pos(2)+21; % something strange here
    a = a(size(a,1)-pos(2)-pos(4)+(0:pos(4)-1),pos(1)+(0:pos(3)-1),:);
    imwrite(a,fname)
    
else
    
    % Linux: First get figure name
    figure(hf), pause(.5)
    if strcmp(get(hf,'numbertitle'),'on')
        if mod(hf,1)
            figname = sprintf('Figure %.6f',hf);
        else
            figname = sprintf('Figure %i',hf);
        end
    else
        figname = '';
    end
    name = get(hf,'name');
    if ~isempty(name)
        if isempty(figname)
            figname = name;
        else
            figname = [figname ': ' name];
        end
    end
    
    % Get figure id using linux command 'xwininfo'
    [status result] = system(['xwininfo -name "' figname '"']);
    if status, disp(result), error('an error occured'), end
    token = regexp(result,'Window id: 0x([a-h\d])*','tokens');
    id = token{1}{1};
    
    % Create file name
    if nargin<2
        fname = [fn_cd('capture') '/' figname  datestr(now,'_yymmdd_HH:MM:SS') '.png'];
    end
    
    % Save figure using ImageMagick function 'import' through the id
    [status result] = system(['import -window 0x' id ' "' fname '"']);
    if status, disp(result), error('an error occured'), end
    
end
