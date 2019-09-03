function fn_progress(varargin)
% function fn_progress(prompt,max[,ht])
% function fn_progress(prompt,max,'%' or 'p'[,ht])
% function fn_progress(prompt)
% function fn_progress(i)
% function fn_progress(i,'pause')
% function fn_progress('end')
% function fn_progress('cont')
% function fn_progress('in',ht)
% function fn_progress('screen')
% function fn_progress('elapsed|elapsedmin')
%---
% progress indicator

% Thomas Deneux
% Copyright 2003-2012

persistent x    % structure with persistent information
persistent ht0  % default place for displaying progress (handle or '' for command prompt)
if ischar(varargin{1}) 
    prompt = varargin{1};
    
    % SPECIAL COMMANDS
    switch prompt
        case 'end'
            if isempty(x), return, end
            if ishandle(x(1).ht)
                set(x(1).ht,'string','')
            else
                fprintf(repmat('\b',1,x(1).promptsize+1+x(1).isize+length(x(1).after)+1))
            end
            x(1) = [];
            return
        case 'cont'
            fprintf(repmat(' ',1,x(1).promptsize+1+x(1).isize+length(x(1).after)))
            return
        case 'elapsed'
            disp(['elapsed ' num2str(toc(x(1).timer)) 's.'])
            return
        case 'elapsedmin'
            disp(['elapsed ' num2str(toc(x(1).timer)/60,'%.1f') 's.'])
            return
        case 'in'
            ht0 = varargin{2};
            return
        case 'screen'
            if ~isempty(x) && any(ishandle(x(1).ht)), set(x(1).ht,'string',''), end
            ht0 = [];
            return
    end
    
    % INITIALIZATION
    % set current structure and remember last parameters in case of nested calls to fn_progress
    x0 = struct( ...
        'prompt',[], ...       % prompt 
        'promptsize',[], ...   % size of prompt - used for 'end' and 'cont'
        'curi',[], ...         % remembers current progress - avoids printing again the same thing
        'after',[], ...        % message after the progress - e.g. %, or /139
        'isize',[], ...        % length of the string for progress number
        'format',[], ...       % format for the number
        'doerase',[], ...      % do erase previous progress - on by default, at init use negative 'max' number to set it off)
        'pflag',[], ...        % use % rather than /139
        'max',[], ...          % maximal progress number
        'ht',[], ...           % current place for displaying progress - if specified at init, is used instead of ht0
        'timer',[], ...        % timer started at init - used for 'elapsed' and 'elapsedmin'
        'caller',[] ...        % name of caller function ('' for base workspace) - used to prevent nested calls to fn_progress
        );
    if isempty(x), x = x0; else x = [x0 x(1)]; end
    % use caller to detect nested calls
    stack = dbstack;
    if isscalar(stack), x(1).caller = ''; else x(1).caller = stack(2).name; end
    % sizes
    x(1).promptsize = length(prompt);
    i = 0; x(1).curi = i;
    if nargin<2
        % use default parameters
        x(1).pflag = true;
        x(1).max = 1;
        x(1).isize = 3;
    else
        x(1).max = varargin{2};
        if ischar(x(1).max)
            x(1).max = str2double(x(1).max); 
        end
        x(1).pflag = false;
        x(1).max = abs(x(1).max);
        x(1).isize = floor(log10(x(1).max))+1;
        if nargin>2
            x(1).ht = varargin{3};
        else
            x(1).ht = ht0;
        end
        x(1).doerase = isempty(x(1).ht) && (x(1).max>0);
    end
    % format
    x(1).prompt = prompt;
    x(1).format = ['%' num2str(x(1).isize) 'i'];
    if x(1).pflag
        x(1).after = sprintf('%%');
    else
        x(1).after = ['/' sprintf(x(1).format,x(1).max)];
    end
    % initial display
    if x(1).doerase, updatedisplay(-1,[prompt ' ' sprintf(x(1).format,0) x(1).after],x(1).ht), end
    % start x(1).timer
    x(1).timer = tic;
else
    % STATE PROGRESS
    i = varargin{1};
    % detect nested call to fn_progress
    stack = dbstack;
    if isscalar(stack), caller = ''; else caller = stack(2).name; end
    if ~strcmp(caller,x(1).caller)
        % double nesting
        if isscalar(x) || ~strcmp(caller,x(2).caller)
            if isscalar(x), disp('fn_progress: strange! unregistered caller'), end
            % here we lost track of the right parameters
            str = ['???? ' num2str(i) '/??'];
            updatedisplay(-1,str,'')
            return
        end
        % cancel the most nested call!
        x(1) = [];
        x(1).doerase = false;
    end
    % display
    if x(1).pflag, i = floor(i/x(1).max*100); end
    if (i == x(1).curi), return, end
    x(1).curi = i;
    if x(1).doerase
        nerase = x(1).isize+length(x(1).after);
        updatedisplay(nerase,[sprintf(x(1).format,i) x(1).after],x(1).ht)
    else
        updatedisplay(-1,[x(1).prompt ' ' sprintf(x(1).format,i) x(1).after],x(1).ht)
    end
    if nargin>1 && strcmp(varargin{2},'pause')
        pause
    end
end

drawnow

%---
function updatedisplay(nerase,str,ht)

if ~isempty(ht) && ishandle(ht)
    if nerase>0
        s = get(ht,'string');
        ns = length(s);
        s = s(1:ns-nerase);
    else
        s = [];
    end
    set(ht,'string',[s str])
    drawnow
else
    fprintf([repmat('\b',1,nerase+1) str '\n'])
end
