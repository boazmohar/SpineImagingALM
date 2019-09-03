function filename = fn_getfile(varargin)
% function filename = fn_getfile(['READ|SAVE',][filter[,title]])
% function filename = fn_getfile('DIR',title)
% function rep = fn_getfile('REP')
% function fn_getfile('REP',rep)
%---
% returns file name with full path
% fn_getfile('REP') returns the current directory
% fn_getfile('REP',rep) sets the current directory
% as a utility, fn_getfile('cd') is going to the current directory 
% 
% See also fn_savefile

% Thomas Deneux
% Copyright 2003-2012

persistent rep

% cd to persistent directory
swd=pwd;
if isempty(rep), rep = pwd; end
try cd(rep), catch, rep = pwd; end

% Input
arg = varargin;
if (length(arg)>=1 && ischar(arg{1}) && fn_ismemberstr(arg{1},{'READ','SAVE','DIR','REP','cd'}))
    flag=arg{1}; 
    arg={arg{2:end}}; 
else
    flag='READ';
end
arg0 = arg;
if length(arg)<1 || isempty(arg{1}), arg{1}='*.*'; end
if length(arg)<2 || isempty(arg{2}), arg{2}='Select file'; end

% different actions
switch flag
    case 'READ'
        [filename pathname] = uigetfile(arg{:},'MultiSelect','on');
        if iscell(filename), filename = strvcat(filename{:}); end %#ok<VCAT>
        if filename
            filename = [repmat(pathname,size(filename,1),1) filename];
            rep = pathname;
        end
    case 'SAVE'
        [filename pathname] = uiputfile(arg{:});
        if filename
            filename = fullfile(pathname,filename);
            rep = pathname;
        end
    case 'DIR'
        if length(arg0)>1, error('too many arguments'), end
        arg = [{[]} arg0];
        if length(arg)<2 || isempty(arg{2}), arg{2}='Select directory'; end
        filename = uigetdir(arg{:});
        if filename
            rep = filename;
        end
    case 'REP'        
        % special case: access/set the persistent directory name rep
        if nargin==1
            filename = rep;
        elseif nargin==2
            rep = varargin{2};
        end
    case 'cd'
end

% cd back to initial directory
if ~strcmpi(flag,'cd'), cd(swd), end