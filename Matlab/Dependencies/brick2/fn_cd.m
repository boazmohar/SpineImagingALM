function repout = fn_cd(flag,varargin)
% function rep = fn_cd(flag,relpath)
% function fn_cd('edit')
%---
% This function allows user to rapidly change the current directory or get
% the full path to some files. 
%
% - Type 'fn_cd edit' to launch a GUI interface that enables you to define
%   your own absolute and relative paths (for example the full path to your
%   home, to which will be associated the flag 'home')
% - Then type 'fn_cd home' to set Matlab current directory to your home.
% - Type 'a = fn_cd('home');' to get the full path to your home (note that
%   this will not change Matlab current directory
% - Additional arguments to fn_cd can be used to access subdirectories or
%   files inside a flagged folder. For example, type 'fn_cd home dir1 dir2'
%   to set Matlab current directory to the subdirectory dir1/dir2 inside
%   your home directory, or type 'a = fn_cd('home','dir1/myfile');' to get
%   the full path to the file 'myfile' inside the subdirectory 'dir1'
%   inside your home directory.

% Thomas Deneux
% Copyright 2002-2012

if nargin==0, help fn_cd, return, end

% User define of new flags
if strcmp(flag,'edit')
    fncdgui
    return
end
    
s = loaddef;
rep = getdir(s,flag);

if nargin>1
    rep = fullfile(rep,varargin{:});
end

if nargout==0
    if ~isempty(rep), cd(rep), end
else 
    repout = rep;
end

end

%---
function hostname = gethostname()
hostname = getenv('HOSTNAME');
if isempty(hostname)
    [dum hostname] = system('echo $HOSTNAME');
    hostname = strrep(hostname,char(10),''); % remove endlines
end %#ok<*ASGLU>
hostname = [computer '-' hostname];
end

%---
function path = resolvehost(path)
if isstruct(path)
    khost = strcmp(gethostname(),{path.host});
    if ~isscalar(khost), path = ''; return, end % error!
    path = path.path;
end
end

%---
function fname = getdir(s,i)
if ischar(i)
    label = i;
    i = find(strcmp(label,{s.label}));
    if ~isscalar(i), fname = ''; return, end % error!
end
if isempty(s(i).relto)
    base = '';
else
    base = getdir(s,s(i).relto);
end
path = resolvehost(s(i).path);
fname = fullfile(base,path);
end

%---
function s = loaddef
fname = [which('fn_cd') 'at'];
if exist(fname,'file')
    load(fname)
else
    s = struct('label',{'example1' 'example2'},'relto',{'' 'example1'},'path',{'/home/deneux' 'Matlab'});
end
end

%---
function savedef(s) %#ok<INUSD>
fname = [which('fn_cd') 'at'];
save(fname,'s')
end

%---
function fncdgui

% position parameters
H = 600;
W = 560;
d = 5;

w1 = 80;
w2 = w1;
w3 = 18;
w3label = 120;
w5 = 35;
w4 = W - (d+w1+d+w2+d+w3+d+d+w5+d);
wp = 150;
ws = 60;
hh = 18;
ht = 23;
hd = ht-hh;
nlin = floor((H-hd)/ht);
ndefperpage = nlin-1;

% figure
hf = figure(726);
fn_setfigsize(hf,W,H);
set(hf,'numbertitle','off','name',['FN_CD (current host: ' gethostname() ')'])
set(hf,'defaultuicontrolhorizontalalignment','left')

% load definitions if exist
s = loaddef();
ndef = length(s);
npage = ceil(ndef/ndefperpage);

% first line
uicontrol('style','text','string','label','pos',[d H-ht w1 hh])
uicontrol('style','text','string','relative to','pos',[d+w1+d H-ht w2 hh])
uicontrol('style','text','string','host dependent?','pos',[d+w1+d+w2+d H-ht w3label hh])
pagestr = [cellstr(num2str((1:npage)','page %i')); 'new page'];
kpage = npage;
hupage = uicontrol('style','popupmenu', ...
    'string',pagestr,'value',kpage, ...
    'pos',[W-d-ws-d-wp H-ht wp hh], ...
    'callback',@(u,e)chgpage(get(u,'value')));
uicontrol('style','pushbutton','string','SAVE', ...
    'pos',[W-d-ws H-ht ws hh], ...
    'callback',@(u,e)saveit())

% definition lines
hu = zeros(3,ndefperpage);
hm = zeros(1,ndefperpage);
for k=1:ndefperpage
    % we first display all the 'edit' controls, so that they come next to
    % each other when pressing the Tab key
    hu(1,k) = uicontrol('style','edit', ...
        'pos',[d H-(k+1)*ht w1 hh], ...
        'callback',@(u,e)chgdef(k));
    hu(2,k) = uicontrol('style','edit', ...
        'pos',[d+w1+d H-(k+1)*ht w2 hh], ...
        'callback',@(u,e)chgdef(k));
    hu(3,k) = uicontrol('style','edit', ...
        'pos',[d+w1+d+w2+d+w3+d H-(k+1)*ht w4 hh], ...
        'callback',@(u,e)chgdef(k));
end
for k=1:ndefperpage
    % we next display the checkboxes (for hostname depedency) and the
    % buttons for user access to directories
    hm(k) = uicontrol('style','checkbox', ...
        'pos',[d+w1+d+w2+d H-(k+1)*ht w3 hh], ...
        'callback',@(u,e)chgdef(k));
    uicontrol('style','pushbutton','string','D', ...
        'pos',[d+w1+d+w2+d+w3+d+w4+d H-(k+1)*ht w5 hh], ...
        'callback',@(u,e)userdir(k));
end
displaypage()

    function displaypage()
        for k=1:ndefperpage
            kdef = (kpage-1)*ndefperpage + k;
            if kdef>ndef || isempty(s(kdef).label)
                % no definition exists
                set(hu(:,k),'string','','backgroundcolor','default')
                set(hm(k),'value',0)
            else
                % definition exists
                sk = s(kdef);
                fn_set(hu(:,k),'string',{sk.label,sk.relto,resolvehost(sk.path)})
                set(hm(k),'value',isstruct(sk.path))
                checkdir(k)
            end
        end
    end
    function chgpage(k)
        kpage = k;
        if kpage>npage
            npage = kpage;
            pagestr = [cellstr(num2str((1:npage)','page %i')); 'new page'];
            set(hupage,'string',pagestr);
        end
        displaypage();
    end
    function chgdef(k)
        kdef = (kpage-1)*ndefperpage + k;
        c = fn_get(hu(:,k),'string');
        [s(kdef).label s(kdef).relto path] = deal(c{:});
        if get(hm(k),'value')
            % host name - specific definition
            sk = s(kdef);
            path = struct('host',gethostname(),'path',path);
            if isstruct(sk.path)
                khost = strcmp(path.host,{sk.path.host});
                if isempty(khost)
                    khost = length(sk.path)+1;
                elseif ~isscalar(khost)
                    error strange
                end
                s(kdef).path(khost) = path;
            else
                s(kdef).path = path;
            end
        else
            % not host name - specific
            s(kdef).path = path;
        end
        checkdir(k)
        if kdef>ndef, ndef=kdef; end
    end
    function checkdir(k)
        kdef = (kpage-1)*ndefperpage + k;
        fname = getdir(s,kdef);
        if isempty(s(kdef).label)
            set(hu(:,k),'backgroundcolor','default','string','')
        elseif exist(fname,'dir')
            set(hu(:,k),'backgroundcolor','w')
        else
            set(hu(:,k),'backgroundcolor',[1 .5 0])
        end
    end
    function userdir(k)
        kdef = (kpage-1)*ndefperpage + k;
        path = fn_getdir;
        if kdef<=ndef && ~isempty(s(kdef).relto)
            s(kdef).path = '';
            relto = getdir(s,kdef);
            if length(path)<length(relto) || ~all(relto==path(1:length(relto)))
                errordlg('selected directory is not inside ''%s''',s(kdef).relto)
                return
            end
            path = path(length(relto)+1:end);
        end
        set(hu(3,k),'string',path)
        chgdef(k)
    end
    function saveit()
        savedef(s)
        assignin('base','s',s)
    end

end










