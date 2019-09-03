function createhelpfiles(mfile)
% function createhelpfiles(mfile)

if nargin==0, help createhelpfiles, return, end

fn_cd brick
if strcmp(mfile,'all')
    d = cellstr(fn_ls('*.m'));
else
    d = {[mfile '.m']};
end
opt = struct('outputDir',fn_cd('brick','html'));

for kfile=1:length(d)
    
    % file names
    mfile = d{kfile};
    if strcmp(mfile,'Contents.m'), continue, end
    disp(mfile)
    if strcmp(mfile,'Contents.m'), continue, end
    hfile = ['helpsource/autohelp/help_' mfile];
    
    % read file
    content = cellstr(fn_readtext(mfile));
    
    % get desc and copyright blocks
    khelpend = 0;
    syntax = {};
    desc = {};
    copyright = {};
    for k=2:length(content)
        line = content{k};
        istart = find(line~=' ',1,'first');
        if isempty(istart) || line(istart)~='%'
            if k==2
                error 'no desc found'
            elseif ~khelpend
                khelpend = k;
            elseif k==khelpend+1
                error 'no copyright found'
            else
                break
            end
            continue
        else
            line = line(istart:end);
        end
        if ~khelpend && isempty(desc) && ~isempty(strfind(line,'function'))
            syntax{end+1} = strrep(line,'function','');
        elseif ~khelpend
            if length(line)<2 || line(2)~='-' % avoid the '%---' line
                desc{end+1} = line;
            end
        else
            copyright{end+1} = line;
        end
    end
    
    % formatting syntax - no need for more formatting
    
    % formatting description
    for k=1:length(desc)
        line = desc{k};
        desc{k} = ['% ' line(2:end)]; % syntax for formatted block
    end
    
    % whole formatted help
    out = {['%% ' fn_fileparts(mfile,'base')]};
    if ~isempty(syntax)
        out = [out ...
            {''} ...
            '%% Syntax' ...
            syntax]; %#ok<*AGROW>
    end
    if ~isempty(desc)
        for k=1:length(syntax)
            line = syntax{k};
            istart = 1+find(line(2:end)~=' ',1,'first');
            syntax{k} = ['%  ' line(istart:end)]; %#ok<*SAGROW> % syntax for formatted block
        end
        out = [out ...
            {''} ...
            '%% Description' ...
            desc];
    end
    if ~isempty(copyright)
        [copyright{2,:}] = deal('%');
        out = [out ...
            {''} ...
            '%% Source' ...
            copyright(:)'];
    end
        
    % finish
    fn_savetext(out,hfile)
    publish(hfile,opt);
    
end