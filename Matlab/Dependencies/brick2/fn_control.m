classdef fn_control < hgsetget
    % function X = fn_control(s[,fun][,spec][,hparent][,'okbutton|nobutton')
    % function fn_control('demo')
    %---
    % cluster of controls representing the structure s (changing the
    % control changes the values in X, and changing the values in X
    % automatically updates the control display)
    %
    % Input:
    % - s           structure to intialize X
    % - fun         function with prototype @(s)fun, which will be called
    %               by X
    % - spec        structure with additional information on the aspect and
    %               behavior of the controls (see below); it should have
    %               the same fields as s (though some fields can be
    %               omitted)
    % - hp          parent figure or uipanel where to create the controls
    %               (a new figure is created if not specified)
    % - 'okbutton' or 'nobutton'
    %               specifically specify to have an ok button or no button
    %
    % Output:
    % - X           a fn_control object, which can be manipulated using
    %               usual structure syntax
    % 
    % Possible values for the fields of spec:
    % []            automatic guess how to display the control
    % 'logical'     check box
    % 'multcheck [n]' array of n check boxes
    % {'str1','str2',...}       
    %               list box with specified string values
    % {'list|radio|button' 'str1','str2',...}       
    %               specification of the type of display [default: list]
    %               for a choice between string values
    % 'char [n [nlin]]'    
    %               input for string, if n is specified, sets the minimal
    %               length of the input in number of characters, otherwise,
    %               minimal length is set according to the value in s
    %               if nlin is specified, control occupies nlin lines
    %               instead of 1
    % 'double [n]'  input for numerical array
    % 'single [n]'  input for numerical array
    % 'slider min max [step] [format]'
    %               slider, specify min, max, step (optional) and format of the
    %               string representation (optional)
    % 'logslider min max [step] [format]'
    %               logarithmic scale slider (min and max should be the log of
    %               the effective min and max)
    % 'loglogslider min max [step] [format]'
    %               logarithmic scale slider, with possibility to select
    %               also a negative number
    % 'stepper [n [min [max [step [format]]]]]'
    %               input for n double
    %               if n>1, it is possible to define n values for min, max,
    %               step, separated by commas, for example: 0,-Inf,-1
    % 'clip'        input for 2-elements vector (usually, min and max);
    %               move the mouse in the control area to change the value
    % 'xdouble, xsingle, xchar [n], x[log[log]]slider min max [..],
    % xstepper, xclip'
    %               additional display of a check box: value will be empty
    %               if the box is not checked
    %               it is possible to specify a default value inside
    %               brackets at the end of the flag, for example:
    %               'xchar 12 [yes we can]' (here the brackets do not mean
    %               that this default value is optional, but they must
    %               appear in the string)
    % 'file|dir'    button for selecting file name / directory name
    % 
    % See also fn_structedit
    
    % Thomas Deneux
    % Copyright 2007-2012
    
    properties
        fun
        immediateupdate
    end
    properties (SetAccess='private')
        mode              % which special buttons: 'none', 'ok', 'execfun' 
        controls
        names
        entries
        hp
        himupd
        fignew
        changedfields
    end
    
    properties (Dependent)
        s
    end
    
    events
        OK
    end
    
    % Constructor
    methods
        function X = fn_control(s,varargin)
            
            % demo
            if ischar(s) && strcmp(s,'demo')
                [s varargin] = demo;
            end
            
            % Input
            X.hp = []; spec = struct;
            % (original structure)
            if ~isscalar(s)
                if length(s)==2
                    spec = s(2); s = s(1);
                else
                    error('input structure must be scalar or have 2 elements')
                end
            end
            % (other inputs)
            for i=1:length(varargin)
                a = varargin{i};
                if ischar(a)
                    % flag for the presence of an 'ok' button
                    X.mode = fn_switch(a, ...
                        {'ok' 'okbutton'},  'ok', ...
                        'nobutton',         'none');
                elseif isa(a,'function_handle')
                    % callback function
                    X.fun = a;
                elseif isstruct(a)
                    % specifications
                    spec = a;
                elseif ishandle(a) || fn_isfigurehandle(a)
                    % parent object (figure or panel)
                    X.hp = a;
                else 
                    error argument
                end
            end
            % (some initializations)
            if isempty(X.mode)
                X.mode = fn_switch(isempty(X.fun),'none','execfun');
            end
            X.fignew = isempty(X.hp);
            X.names = fieldnames(s)';
            X.entries = struct;
            nf = length(X.names);
            
            % Make a general structure with a lot of information
            X.controls = struct('name',{},'type',{},'style',{},'value',{}, ...
                'default',{},'defaultcheck',{},'defaultstring',{}, ...
                'n_name',{},'n_val',{},'n_line',{},'hname',{},'hval',{}, ...
                'check',{},'log',{},'min',{},'max',{},'step',{},'shift',{},'format',{}, ...
                'mode',{},'values',{} ...
                );
            
            % Loop on fields: fill the structure with all necessary
            % parameters but do not display the controls yet
            for k=1:nf
                f = X.names{k};
                X.entries.(f) = k;
                X.controls(k).name = f;
                xk = X.controls(k);
                xk.value = s.(f);
                if isfield(spec,f), opt = spec.(f); else opt = []; end

                % field name
                xk.name = f;
                xk.n_name = length(f);
                
                % default specification when not specified
                if isempty(opt)
                    opt = fn_switch(class(xk.value), ...
                        'logical',  'logical', ...
                        'double',   'double', ...
                        'single',   'single', ...
                        'char',     'char');
                end
                
                % value type and control style
                if iscell(opt)
                    % (list of strings)
                    xk.type = 'char';
                    xk.check = false;
                    if ismember(opt{1},{'list' 'radio' 'button'})
                        % format {'style','str1','str2',...}
                        xk.style = fn_switch(opt{1}, ...
                            'list',     'popupmenu', ...
                            'radio',    'radiobutton', ...
                            'button',   'togglebutton');
                        opt(1) = []; spec.(f) = opt;
                    else
                        % format {'str1','str2',...}, default style is 'list'
                        xk.style = 'popupmenu';
                    end
                    xk.defaultcheck = true;
                else
                    % (string describing the control)
                    % check first flag
                    tmp = regexp(opt,'^[^ ]*','match'); tmp = tmp{1};
                    % check box?
                    if tmp(1) == 'x'
                        xk.check = true;
                        tmp(1) = [];
                        % initial checking of the box
                        xk.defaultcheck = ~isempty(xk.value);
                    else
                        xk.check = false;
                        xk.defaultcheck = true;
                    end
                    % define value type and display style
                    switch tmp
                        case 'logical'
                            xk.type = 'logical';
                            xk.check = true;
                            xk.style = 'checkbox';
                        case 'multcheck'
                            xk.type = 'logical';
                            xk.style = 'multcheck';
                        case {'slider' 'logslider' 'loglogslider'}
                            xk.type = 'double';
                            xk.log = sum(logical(findstr(tmp,'log'))); % count how many 'log'
                            xk.style = 'slider';
                        case 'stepper'
                            xk.type = 'double';
                            xk.style = 'stepper';
                        case 'clip'
                            xk.type = 'double';
                            xk.style = 'sensor';
                            xk.mode = 'clip';
                        case {'file','dir'}
                            xk.type = 'char';
                            xk.style = 'file';
                            xk.mode = fn_switch(tmp,'file','READ','dir','DIR');
                        case {'char' 'double' 'single'}
                            xk.type = tmp;
                            xk.style = 'edit';
                        otherwise
                            error('unknown control flag: ''%s''',tmp)
                    end
                    % initial (default) value for unchecked check box
                    idx = find(opt=='[');
                    if ~xk.defaultcheck && ~isempty(idx)
                        % for the moment, we keep a string for the
                        % default value, which will be converted later
                        % to the appropriate type
                        str = opt(idx+1:end-1);
                        xk.default = str2val(str,xk.type);
                    end
                    opt(idx:end) = [];
                end
                
                % initial value, control width, and some other
                % style-specific parameters
                xk.n_line = 1;
                switch xk.style
                    case 'checkbox'
                        % logical value - only check box
                        xk.n_val = 0;
                        if isempty(xk.value)
                            xk.default = false;
                        else
                            xk.default = logical(xk.value);
                        end
                    case 'multcheck'
                        % read options
                        answer = regexp(opt,'([^ ]*)','tokens');
                        if isempty(xk.value)
                            if length(answer)==2
                                xk.mode = str2double(answer{2}); % number of check boxes
                            else
                                xk.mode = 1;
                            end
                            xk.default = false(1,xk.mode);
                        else
                            if length(answer)==2 && str2double(answer{2})~=length(xk.value)
                                error('multcheck length specification does not match with value')
                            end
                            xk.mode = length(xk.value);
                            xk.default = xk.value;
                        end
                        % width
                        xk.n_val = 3*xk.mode; % not clear what is the minimum width
                        xk.format = '%i';
                    case {'popupmenu' 'radiobutton' 'togglebutton'}
                        % (list of strings)
                        if isempty(xk.value)
                            xk.default = 1;
                        else
                            tmp = find(strcmp(xk.value,opt));
                            if isempty(tmp)
                                xk.default = 1;
                            else
                                xk.default = tmp;
                            end
                        end
                        xk.values = opt;
                        tmp = strvcat(opt{:}); %#ok<VCAT>
                        switch xk.style
                            case 'popupmenu'
                                xk.n_val = size(tmp,2)*.7;
                            case 'togglebutton'
                                xk.n_val = length(opt) + numel(tmp)*.7;
                            case 'radiobutton'
                                xk.n_val = 4*length(opt) + numel(tmp)*.7;
                        end
                    case 'slider'
                        % read options
                        answer = regexp(opt,'([^ ]*)','tokens');
                        answer = [answer{:}];
                        if length(answer)>=3
                            xk.min = str2double(answer{2});
                            xk.max = str2double(answer{3});
                        else
                            xk.min = fn_switch(xk.value>0,0,2*xk.value);
                            xk.max = fn_switch(xk.value==0,1,2*abs(xk.value));
                        end
                        xk.step = 0;
                        xk.format = [];
                        for i=4:length(answer)
                            str = answer{i};
                            if findstr(str,'%')
                                xk.format = str;
                            else
                                xk.step = str2double(str);
                            end
                        end
                        if isempty(xk.format)
                            if xk.step>0 && ~mod(xk.step,1)
                                xk.format = '%.0f';
                            else
                                xk.format = '%.1f';
                            end
                        end
                        if xk.log==2
                            xk.shift = xk.min;
                            xk.max = xk.max-xk.min;
                            xk.min = -xk.max;
                        end
                        % initial value
                        if xk.defaultcheck
                            if ~isnumeric(xk.value) || ~isscalar(xk.value)
                                error('wrong default value in the context of  slider style')
                            end
                            switch xk.log
                                case 0
                                    xk.default = xk.value;
                                case 1
                                    xk.default = log10(xk.value);
                                case 2
                                    % complicate
                                    xk.default = x2log(xk.value,xk.shift);
                            end
                        elseif isempty(xk.default)
                            if xk.log==2
                                xk.default = 0;
                            else
                                xk.default = xk.min;
                            end
                        end
                        % increase the field length according to format
                        test1 = num2str(xk.max+rand,xk.format);
                        test2 = num2str(xk.min+rand,xk.format);
                        xk.n_name = xk.n_name + 2 + max(length(test1),length(test2)) + 1;
                        % second column
                        xk.n_val = 6;
                    case 'stepper'
                        % read options
                        defans = {'stepper' '1' '-Inf' 'Inf' '1' '%.2g'};
                        answer = regexp(opt,'([^ ]*)','tokens');
                        answer = [answer{:}];
                        missing = (length(answer)+1:length(defans));
                        answer(missing) = defans(missing);
                        xk.mode = str2double(answer{2}); % number of numeric values
                        xk.min  = str2num(['[' strrep(answer{3},',',' ') ']']); %#ok<ST2NM>
                        xk.max  = str2num(['[' strrep(answer{4},',',' ') ']']); %#ok<ST2NM>
                        xk.step = str2num(['[' strrep(answer{5},',',' ') ']']); %#ok<ST2NM>
                        xk.format = answer{6};
                        % initial value
                        if xk.defaultcheck
                            xk.default = xk.value; 
                        elseif isempty(xk.default)
                            xk.default = zeros(1,xk.mode);
                        end
                        % width
                        xk.n_val = 5*xk.mode; % not clear what is the minimum width
                    case 'sensor'
                        if xk.defaultcheck
                            xk.default = xk.value; 
                        elseif isempty(xk.default)
                            xk.default = [0 1]; % TODO: check this
                        end
                        xk.n_val = 10; % not clear what is the minimum width
                    case 'edit'
                        if xk.defaultcheck
                            xk.default = val2str(xk.value); 
                        elseif isempty(xk.default)
                            xk.default = '';
                        end
                        tmp = regexp(opt,'[^ ]*','match');
                        if length(tmp)>=2 % length is specified
                            xk.n_val = str2double(tmp{2});
                        else
                            xk.n_val = max(4,length(xk.default));
                        end
                        if strcmp(xk.type,'char') && length(tmp)>=3
                            nlin = str2double(tmp{3});
                            if nlin<=0 || mod(nlin,1), error('number of lines must be a positive integer'), end
                            xk.n_line = nlin;
                            xk.n_val = xk.n_val+8;
                        end
                    case 'file'
                        if ~ischar(xk.value)
                            xk.value = '';
                        end
                        xk.default = xk.value;
                        xk.n_val = 20;
                end
                X.controls(k) = xk;
            end
            
            % Width of the two columns
            idx = logical([X.controls.n_val]); % ignore fields which do not have two columns
            if ~X.fignew, fsz = get(X.hp,'defaultuicontrolfontsize'); else fsz=10; end
            if fsz==8 || (fsz==10 && strcmp(computer,'MACI64'))
                n1 = [X.controls(idx).check]*20 + 5 + [X.controls(idx).n_name]*5;
                n2 = 20 + [X.controls(idx).n_val]*3;
            elseif fsz==10
                n1 = [X.controls(idx).check]*20 + [X.controls(idx).n_name]*8;
                n2 = 20 + [X.controls(idx).n_val]*8;
            end
            A = max(n1);
            B = max(n2);
            
            % Position parameters 
            D = 5; E = 10;
            G = 45; L = 80;
            K = 20; T = 17;
            nbut = sum([X.controls.n_line]) + ~strcmp(X.mode,'none');
            if X.fignew
                % figure
                ss = get(0,'screensize');
                W = max(D+A+D+B+D,D+G+D+L+D);
                H = E+nbut*(K+E);
                H = H+20; % BUG with Exceed
                p = get(0,'pointerLocation');
                LEFT   = min(max(p(1)-W/2,20),ss(3)-W-20);
                BOTTOM = min(max(p(2)-H/2,50),ss(4)-H-20);
                posp = [LEFT BOTTOM W H];
                X.hp = figure('numbertitle','off','name','Edit parameters', ...
                    'menubar','none', ...
                    'position',posp, ...
                    'defaultuicontrolhorizontalalignment','left');
                ncol = 1;
                nlin = nbut;
            else
                % check size of parent
                switch get(X.hp,'type')
                    case 'figure'
                        posp = get(X.hp,'position');
                    case 'uipanel'
                        oldunit = get(X.hp,'units');
                        set(X.hp,'units','pixel')
                        posp = get(X.hp,'position');
                        set(X.hp,'units',oldunit);
                    otherwise
                        error('containter must be a figure or a uipanel object')
                end
                Z = max(A+D+B+D, ...                                        % maximal width of normal buttons
                    fn_switch(X.mode,'execfun',G+D+L+D,'ok',G+D,'none',0)); % width of special button
                Y = K+E;
                xrep = (posp(3)-D)/Z;
                yrep = (posp(4)-D)/Y;
                if floor(xrep)*floor(yrep)>=nbut
                    % it fits
                    ncol = ceil(nbut/floor(yrep));
                    nlin = ceil(nbut/ncol);
                    B = (posp(3)-D-ncol*(A+D+D))/ncol;
                elseif floor(xrep)*floor(yrep/.7)>=nbut
                    % it fits when squeezed a bit vertically
                    ncol = floor(xrep);
                    nlin = ceil(nbut/ncol);
                    B = (posp(3)-D-ncol*(A+D+D))/ncol;
                    K = 17;
                    E = (posp(4)-nlin*K)/(nlin+1);
                else
                    error('cannot fit the buttons inside the container')
                end
            end
%             delete(get(X.hp,'children'))
%             set(X.hp,'DeleteFcn',@(hp,evnt)delete(X))
            
            bgcol = get(X.hp,fn_switch(get(X.hp,'type'),'figure','color','uipanel','backgroundcolor'));
            set(X.hp,'defaultuicontrolunits','pixel')
            X.changedfields = {};
            H = posp(4);
            Z = max(A+D+B+D, ...                                        % maximal width of normal buttons
                fn_switch(X.mode,'execfun',G+D+L+D,'ok',G+D,'none',0)); % width of special button
            Y = K+E;
            
            % Controls
            ipos = 1;
            for k=1:nf
                xk = X.controls(k);
                icol = 1+floor((ipos-1)/nlin);
                ilin = ipos-(icol-1)*nlin;
                ipos = ipos + xk.n_line; % position for next control
                if ~xk.n_val
                    % first column only (checkbox)
                    xk.hname = uicontrol('style','checkbox','string',xk.name, ...
                        'value',xk.default, ...
                        'parent',X.hp,'backgroundcolor',bgcol, ...
                        'position',[D+(icol-1)*Z H-ilin*Y A+D+B T], ...
                        'callback',@(hu,evnt)chgvalue(X,k,logical(get(hu,'value'))));
                else
                    % first column
                    xk.hname = uicontrol('style','text','string',xk.name, ...
                        'horizontalalignment','left', ...
                        'parent',X.hp,'backgroundcolor',bgcol, ...
                        'position',[D+(icol-1)*Z H-ilin*Y A T]);
                    if xk.check
                        set(xk.hname,'style','checkbox','value',xk.defaultcheck, ...
                            'callback',@(hu,evnt)chgvalue(X,k,logical(get(hu,'value'))));
                    end
                    if strcmp(xk.style,'slider') && (~xk.check || xk.defaultcheck)
                        set(xk.hname,'string',[xk.name ' (' num2str(xk.value,xk.format), ')'])
                    end
                    
                    % second column
                    switch xk.style
                        case 'popupmenu'
                            xk.hval = uicontrol('parent',X.hp,'style',xk.style, ...
                                'position',[D+(icol-1)*Z+A+D H-ilin*Y B K], ...
                                'string',xk.values,'value',xk.default, ...
                                'callback',@(hu,evnt)chgvalue(X,k,xk.values{get(hu,'value')}));
                        case {'radiobutton' 'togglebutton'}
                            compactstyle = fn_switch(xk.style, ...
                                'radiobutton',  'radio', ...
                                'togglebutton', 'toggle');
                            xk.hval = fn_buttongroup(compactstyle,xk.values, ...
                                @(x)chgvalue(X,k,x), ...
                                'parent',X.hp, ...
                                'units','pixel','position',[D+(icol-1)*Z+A+D H-ilin*Y B K], ...
                                'value',xk.default);
                        case 'multcheck'
                            xk.hval = fn_multcheck('parent',X.hp, ...
                                'position',[D+(icol-1)*Z+A+D H-ilin*Y B K], ...
                                'value',xk.default,'callback', ...
                                @(hu,evnt)chgvalue(X,k,get(hu,'value')));
                        case 'slider'
                            xk.hval = fn_slider('parent',X.hp,'mode','point', ...
                                'position',[D+(icol-1)*Z+A+D H-ilin*Y B K], ...
                                'min',xk.min,'max',xk.max,'width',.1, ...
                                'value',xk.default,'callback', ...
                                @(hu,evnt)chgvalue(X,k,get(hu,'value')));
                            if xk.step
                                % note that 'inc' and 'width are relative values (betw.
                                % 0 and 1, not betw. min and max)
                                inc = xk.step/(xk.max-xk.min);
                                set(xk.hval,'inc',inc,'width',inc)
                            end
                        case 'stepper'
                            xk.hval = fn_stepper('parent',X.hp, ...
                                'position',[D+(icol-1)*Z+A+D H-ilin*Y B K], ...
                                'min',xk.min,'max',xk.max,'step',xk.step, ...
                                'format',xk.format, ...
                                'value',xk.default,'callback', ...
                                @(hu,evnt)chgvalue(X,k,get(hu,'value')));
                        case 'sensor'
                            xk.hval = fn_sensor('parent',X.hp,'mode',xk.mode, ...
                                'backgroundcolor',[.5 .6 .6], ... 
                                'position',[D+(icol-1)*Z+A+D H-ilin*Y B K], ...
                                'value',xk.default, ...
                                'format','%.4g', ...
                                'callback',@(hu,evnt)chgvalue(X,k,get(hu,'value')));
                        case 'edit'
                            xk.hval = uicontrol('parent',X.hp,'style',xk.style, ...
                                'position',[D+(icol-1)*Z+A+D H-(ilin+xk.n_line-1)*Y B K+(xk.n_line-1)*Y], ...
                                'horizontalalignment','left', ...
                                'max',xk.n_line, ... % allow multiple lines if n_line>1
                                'string',val2str(xk.default),'backgroundcolor','w', ...
                                'callback',@(hu,evnt)chgvalue(X,k,str2val(get(hu,'string'),xk.type)));
                        case 'file'
                            xk.hval = uicontrol('parent',X.hp,'style','edit', ...
                                'position',[D+(icol-1)*Z+A+D H-ilin*Y B K], ...
                                'string',xk.default,'enable','inactive', ...
                                'buttondownfcn',@(hu,evnt)chgvalue(X,k,fn_getfile(xk.mode)));
                        otherwise
                            error programming
                    end
                end
                X.controls(k) = xk;
            end
            
            % Special action buttons
            switch X.mode
                case 'execfun'
                    X.immediateupdate = true;
                    uicontrol('style','pushbutton','string','Update','parent',X.hp, ...
                        'position',[D+ncol*Z-D-L-D-G H-nlin*Y G K], ...
                        'callback',@(hu,evnt)evalfun(X))
                    X.himupd = uicontrol('style','radiobutton','string','immediate', ...
                        'value',X.immediateupdate, ...
                        'parent',X.hp,'backgroundcolor',bgcol, ...
                        'position',[D+ncol*Z-D-L H-nlin*Y L K], ...
                        'callback',@(hu,evnt)set(X,'immediateupdate',get(hu,'value')));
                case 'ok'
                    uicontrol('style','pushbutton','string','OK', ...
                        'position',[W/2-G/2 D G K], ...
                        'callback',@(u,evnt)okpress);
                case 'none'
                    % no special buttons
                    X.immediateupdate = ~isempty(X.fun);
            end
            assignin('base','X',X)
            
            % nested function for when the OK button is pressed
            function okpress
                notify(X,'OK')
                % delete the controls now
                if X.fignew
                    close(X.hp)
                else
                    delete(get(X.hp,'children'))
                end
            end
            
        end
    end
    
    % Get/Set
    methods
        function s = get.s(X)
            c = [X.names; {X.controls.value}];
            for k=1:length(X.controls), if iscell(c{2,k}), c{2,k} = {c{2,k}}; end, end %#ok<CCAT1>
            s = struct(c{:});
        end
        function set.immediateupdate(X,value)
            set(X.himupd,'value',value) %#ok<*MCSUP>
            X.immediateupdate = value;
        end        
    end
    
    % Referencing, assignment, display
    methods
        function x = subsref(X,f)
            switch f(1).type
                case '()'
                    x = X(f(1).subs{:});
                case '.'
                    k = strcmp(X.names,f(1).subs);
                    if ~any(k)
                        x = X.(f(1).subs);
                    else
                        x = X.controls(k).value;
                    end
                otherwise
                    error('wrong referencing of fn_control object')
            end
            if length(f)>1, x = subsref(x,f(2:end)); end
        end
        function X = subsasgn(X,f,x)
            switch f(1).type
                case '()'
                    if length(f)>1
                        subsassgn(X(f(1).sub{:}),f(2:end),x)
                    else
                        X(f(1).subs{:}) = x;
                    end
                case '.'
                    k = strcmp(X.names,f(1).subs);
                    if ~any(k)
                        if length(f)>1
                            X.(f(1).subs) = subsassgn(X.(f(1).subs),f(2:end),x);
                        else
                            X.(f(1).subs) = x;
                        end
                    else
                        if length(f)>1
                            X.controls(k).value = subsassgn(X.controls(k).value,f(2:end),x);
                        else
                            X.controls(k).value = x;
                        end
                        updatecontrol(X,k)
                    end
                otherwise
                    error('wrong referencing of fn_control object')
            end
        end
        function disp(X)
            disp(X.s)
        end
    end
    
    % Routines
    methods
        function evalfun(X)
            if ~isempty(X.fun), feval(X.fun,X.s), end
            if isvalid(X) && ishandle(X.hp), X.changedfields = {}; end
        end
        function chgvalue(X,k,val)
            % callback function executed when control k has been moved
            xk = X.controls(k);
            % check the box if we moved the adjacent control
            if xk.check && ~strcmp(xk.style,'checkbox')
                if islogical(val)
                    if val
                        if strcmp(xk.style,'edit')
                            val = str2val(get(xk.hval,'string'),xk.type);
                        else
                            val = get(xk.hval,'value');
                        end
                    else
                        if strcmp(xk.type,'char')
                            val = '';
                        else
                            val = [];
                        end
                    end
                else
                    % check the box
                    set(xk.hname,'value',true)
                end
            end
            % special action
            switch xk.style
                case 'slider'
                    % logarithmic value
                    switch xk.log
                        case 0
                            % value is ok
                        case 1
                            val = 10^val;
                        case 2
                            % complicate
                            val = log2x(val,xk.shift);
                    end
                    % text update
                    if isempty(val)
                        set(xk.hname,'string',xk.name)
                    else
                        set(xk.hname,'string', ...
                            [xk.name ' (' num2str(val,xk.format) ')'])
                    end
                case 'file'
                    if isequal(val,0), return, end
                    set(xk.hval,'string',val)
            end
            % store the information on which fields were changed
            X.changedfields = union(X.changedfields,xk.name);
            % update value and eval function
            X.controls(k).value = val;
            if X.immediateupdate
                evalfun(X)
            end
        end
        function updatecontrol(X,k)
            % update display of control k (normally, upon a change of the
            % related value)
            xk = X.controls(k);
            % check the box according to whether the value is empty
            if xk.check && ~strcmp(xk.style,'checkbox')
                set(xk.hname,'value',~isempty(xk.value))
            end
            % change value display
            switch xk.style
                case 'checkbox'
                    % logical
                    set(xk.hname,'value',xk.value)
                case {'popupmenu' 'radiobutton' 'togglebutton'}
                    val = xk.value;
                    if ischar(val), val = find(strcmp(xk.values,xk.value)); end
                    if isempty(val), error('value ''%s'' does not exist',xk.value), end
                    set(xk.hval,'value',val);
                case 'multcheck'
                    set(xk.hval,'value',xk.value)
                case 'slider'
                    if isempty(xk.value)
                        set(xk.hname,'string',xk.name)
                    else
                        set(xk.hname,'string',[xk.name ...
                            ' (' num2str(xk.value,xk.format), ')'])
                        switch xk.log
                            case 0
                                val = xk.value;
                            case 1
                                val = log10(xk.value);
                            case 2
                                % complicate
                                val = x2log(xk.value,xk.shift);
                        end
                        set(xk.hval,'value',val);
                        % update in case of coercing actions!
                        val = get(xk.hval,'value');
                        switch xk.log
                            case 0
                                X.controls(k).value = val;
                            case 1
                                X.controls(k).value = 10^val;
                            case 2
                                % complicate
                                X.controls(k).value = log2x(val,xk.shift);
                        end
                    end
                case {'stepper','sensor'}
                    if ~isempty(xk.value)
                        set(xk.hval,'value',xk.value)
                    end
                case 'edit'
                    set(xk.hval,'string',val2str(xk.value));
                case 'file'
                    set(xk.hval,'string',xk.value);
                otherwise
                    error programming
            end
        end
    end
        
    % Misc
    methods
        function access(X) %#ok<MANU>
            keyboard
        end
    end
end


%---
function [str type] = val2str(val)

type = class(val);
switch type
    case 'char'
        str = val;
    case {'double','single','logical'}
%         if ndims(val)>2 || numel(val)>15
%             error('cannot display array with more than 2 dimensions or 15 elements')
%         end
        str = num2str(val);
        str(:,end+1) = ';'; str(:,end+1) = ' ';
        str = reshape(str',1,numel(str));
        str(end-1:end) = [];
        str = regexprep(str,' *',' ');
    case 'cell'
        % special! cell array of strings
        type = 'char';
        str = val;
    otherwise
        error('cannot display object of class ''%s''',type)
end

end

%---
function val = str2val(str,type)

switch type
    case 'char'
        val = str;
    case {'double','logical'}
        val = eval(['[' str ']']);
    otherwise
        error programming
end

end

%---
function val = x2log(x,shift)

if x==0
    val=0;
elseif x>0
    val = max(0,log(x)-shift);
else
    val = min(0,-(log(-x)-shift));
end

end

%---
function x = log2x(val,shift)

if val==0
    x = 0;
elseif val>0
    x = 10^(val+shift);
else
    x = -10^(-val+shift);
end

end

%---
function [s args] = demo %#ok<STOUT>

C = {'s = struct(''a'',false,''b'',1,''c'',2,''d'',''hello'',''e'',[0 1],''f'',pwd);'
    'spec = struct(''c'',''xslider 0 10 1'',''d'',{{''hello'',''yo''}},''e'',''clip'',''f'',''dir'');'
    'myfun = @disp;'
    'fn_control(s,spec,myfun);'};
for k=1:4, disp(C{k}), end
for k=1:3, evalin('base',C{k}), eval(C{k}), end
args = {spec myfun};

end