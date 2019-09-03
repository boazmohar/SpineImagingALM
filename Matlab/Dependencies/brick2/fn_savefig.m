function fn_savefig(varargin)
% function fn_savefig([hf][,fnames,...[,methods...]][,scaling])
%---
%
% Input:
% - hf          figure handles
% - fnames      file name where to save the figure
% - methods     for example: 'jpg', 'eps', 'fig', 'png'
% - scaling     scalar: how much to increase the size of the figure (use
%               small values to increase font size)
%
% For example:
% * fn_savefig(1,'myfigure','png','eps') saves figure 1 in files
% 'myfigure.png' and 'myfigure.eps' 
% * fn_savefig(1:3,'a.png','b.jpg','c.fig',2) save figure 1 in 'a.png',
% figure 2 in 'b.jpg' and figure 2 in 'c.fig'; image sizes are double as
% default
% * fn_savefig(1) asks you where to save figure 1

% Thomas Deneux
% Copyright 2003-2012

% Input
nextarg = 1;
% (which figure(s))
if nargin==0
    hfig = findobj('type','figure');
elseif ishandle(varargin{nextarg})
    hfig = varargin{nextarg};
    nextarg = nextarg+1;
else
    hfig = gcf;
end
% (file names and methods)
nfig = length(hfig);
fname = cell(1,nfig);
if nargin>=nextarg && ischar(varargin{nextarg})
    for i=1:nfig
        fname{i} = varargin{nextarg};
        nextarg = nextarg+1;
    end
    k=0; 
    while nargin>=nextarg+k && ischar(varargin{nextarg+k}), k=k+1; end
    methods = varargin(nextarg+(0:k-1));
    nextarg = nextarg+k;
    if isempty(methods)
        [dum1 dum2 ext] = fileparts(fname{1});
        if ~isempty(ext)
            methods = lower(ext(2:end)); 
        else 
            methods = {'eps','png'}; 
        end
    end
else
    for i=1:nfig
        figure(hfig(i))
        fname{i} = input(['sauver ' num2str(hfig(i)) ' dans ? '],'s');
    end
    methods = {'eps','jpg'};
end
% (scaling)
if nargin>=nextarg
    scaling = varargin{nextarg};
    nextarg = nextarg+1;
else
    scaling = 1;
end


% Save
for k=1:nfig
    hf = hfig(k);
    if isempty(fname{k}), continue, end
    % change paper-position property
    %if all(get(hf,'paperposition')==get(0,'defaultfigurepaperposition')) % not changed by user
    pos = get(hf,'position');                       % position in the screen
    set(hf,'paperUnits','inches','paperposition',[1 1 pos([3 4])*(scaling/90)])     % keep the same image ratio
    %end
    [dum1 dum2 ext] = fileparts(fname{k});
    if isempty(ext)
        methk = methods; 
    else
        if isempty(dum1), fname{k} = dum2; else fname{k} = [dum1 filesep dum2]; end
        methk = {ext(2:end)}; 
    end
    for i=1:length(methk)
        switch methk{i}
            case 'eps'
                saveas(hf,[fname{k} '.eps'],'psc2')
            case 'png'
                % do not print UI control objects (they create errors) ->
                % need to call 'print' function with -noui option
                print(hf,[fname{k} '.png'],'-dpng','-noui')
            otherwise
                saveas(hf,[fname{k} '.' methk{i}],methk{i})
        end
    end
end

% reset paperposition property to allow use of saveas


