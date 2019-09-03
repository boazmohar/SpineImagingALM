classdef fn_supercontrol < hgsetget
    % function X = fn_supercontrol([hp,]specs[,callback[,x0]])
    % function specs = fn_supercontrol.examplespecs
    %---
    % Input:
    % - hp          uipanel handle
    % - specs       structure with fields:
    %               . name      string
    %               . controls  structure with fields (* is mandatory):
    %                           style*  popupmenu, edit, checkbox,
    %                                   pushbutton or stepper
    %                           string  the 'string' property of the
    %                                   control
    %                           length* relative width occupied by the
    %                                   control + its label if any
    %                           default* default value (type depends on
    %                                   control style)
    %                           label   a label placed on the left
    %                           labellength     relative width occupied by
    %                                   the label (set to 0 if no label)
    %                           callback        for push button only:
    %                                   function with prototype 
    %                                   valuek = fun(value) where value is
    %                                   a cell array and valuek an element
    %                                   of this array [MORE DOC NEEDED]
    %                           more    more properties stored in a cell
    %                                   array with successive pairs of
    %                                   property names/values
    % - callback    function to be executed when control values are
    %               changed by user, with prototype @(x)fun(x),
    %               where x is X.x (see below)
    % - x0          initial value (see below comments on X.x)
    %
    % Output:
    % - X           fn_supercontrol object; X.x is a structure that
    %               stores the values, with fields:
    %               .name       string
    %               .active     logical
    %               .value      cell array with values (one per
    %                           control in the specs of the same
    %                           name)
    %
    % One can change the values either by user action (acionning
    % the controls) or by setting the value of X.x.
    %
    % Special notes:
    % - chgactive flag: The callback might be invoked but the
    % "active" part of the value has not changed (for example, when
    % a new, inactive line has been creadted). The property
    % X.activechg says whether this active part has changed or not.
    % - 'edit': if a given control has the style 'edit', but its
    % data is set to numeric value instead of a string, the data is
    % stored in the 'userdata' property of the control, and the
    % string 'userdata' is diplayed.
    % - 'pushbutton': when pressing the button, the callback is executed
    % and changes the value accordingly
    %
    % See also fn_control
 
    % Thomas Deneux
    % Copyright 2010-2012
   
    % Content
    properties
        hp
        x_data = struct('name',{},'active',{},'value',{});
        callback
        specs
        activechg = false;  % change in the active part?
    end
    properties (Dependent)
        x
    end
    properties (Dependent, SetAccess='private')
        nx
    end
    % Private 
    properties
        controls = struct( ...
            'panel',{}, ...
            'check',{}, ...
            'contentpanel',{},'content',{}, ...
            'move',{},'close',{});
        addlinecontrol
        ncharname
    end
   
    % Constructor
    methods
        function X = fn_supercontrol(varargin) 
            % function X = fn_supercontrol([hp,]specs[,callback[,x0]])
            % input
            if isstruct(varargin{1})
                hp = gcf;
            else
                hp = varargin{1}; varargin(1) = [];
            end
            specs = varargin{1};
            if length(varargin)<2, callback = []; else callback = varargin{2}; end
            if length(varargin)<3, x0 = []; else x0 = varargin{3}; end
            
            % some default values inside the panel
            set(hp,'defaultuicontrolfontsize',8)
            delete(get(hp,'children'))
            
            % maximum number of character for the names
            set(hp,'units','pixel')
            siz = get(hp,'position'); siz = siz(3);
            allnames = strvcat(specs.name); %#ok<VCAT>
            X.ncharname = size(allnames,2);
            
            % set properties
            X.hp = hp;
            X.callback = callback;
            X.specs = specs;
            
            % default background color
            colprop = fn_switch(get(hp,'type'),'figure','color','uipanel','backgroundcolor');
            defaultcolor = get(hp,colprop);
            try set(hp,'defaultuipanelbackgroundcolor',defaultcolor), end
            set(hp,'defaultuicontrolbackgroundcolor',defaultcolor)
            
            % 'add line' control
            display_addlinecontrol(X)
            
            % default value (automatic display)
            if ~isempty(x0), X.x = x0; end
        end
    end
    
    % Give an example specification
    methods (Static)
        function specs = examplespecs()
            specs = struct( ...
                'name',     'myname', ...
                'controls', struct('style','edit','string','','length',2,'default',{'hip' 'hop'},'label',{'firstname' 'lastname'},'labellength',1,'more',{{'backgroundcolor','y','fontweight','bold'}}) ...
                );
        end
    end
    
    % Get/Set
    methods
        function x = get.x(X)
            x = X.x_data;
        end
        function set.x(X,x)
            oldx = X.x_data;
            if isequal(oldx,x), return, end
            if isempty(x), x = struct('name',{},'active',{},'value',{}); end % just to make sure
            X.x_data = x;
            if length(x)==length(oldx) && isequal({x.name},{oldx.name})
                % change line contents only
                for i=1:length(x), display_linecontent(X,i), end
            else
                % re-display everything
                display_clear(X)
                for i=1:length(x), display_line(X,i), end
            end
        end
        function nx = get.nx(X)
            nx = length(X.x_data);
        end
    end
    
    % Display
    % (all these functions update X.controls and the display; X.x is
    % supposed to be already set to the new value)
    methods
        function [W H h hdec hbut hdectext htext] = display_sizes(X)
            tmp = get(X.hp,'position'); W = tmp(3); H = tmp(4);
            h = 25;
            hdec = 5;
            hbut = h-8;
            hdectext = (h-10)/2+1; 
            htext = 10;
        end
        function display_addlinecontrol(X)
            [W H h hdec hbut] = display_sizes(X);
            X.addlinecontrol = uicontrol('parent',X.hp,'style','popupmenu', ...
                'units','pixel','position',[10 H-h+hdec min(W-10,max(100,X.ncharname*8)) hbut], ...
                'string',{'add line...' X.specs.name}, ...
                'callback',@(u,e)data_addline(X));
        end
        function display_line(X,i)
            [W H h hdec hbut hdectext htext] = display_sizes(X);
            hpi = uipanel('parent',X.hp, ...
                'units','pixel','position',[1 H-(i+1)*h W h], ...
                'userdata',i, ...
                'borderwidth',0);
            X.controls(i).panel = hpi;
            wcheck = 15+7*X.ncharname;
            wplusminus = 2+2*hbut;
            wcontent = W-wcheck-wplusminus;
            kspec = strcmpi(X.x(i).name,{X.specs.name});
            X.controls(i).check = uicontrol('parent',hpi,'style','checkbox', ...
                'units','pixel','position',[1 hdec wcheck hbut], ...
                'string',X.specs(kspec).name, ...
                'value',X.x(i).active, ...
                'callback',@(u,e)data_toggleactive(X,u));
            X.controls(i).contentpanel = uipanel('parent',hpi, ...
                'units','pixel','position',[wcheck+1 1 wcontent h], ...
                'borderwidth',0);
            X.controls(i).move = uicontrol('parent',hpi,'style','slider', ...
                'units','pixel','position',[W-wplusminus+2 hdec hbut-1 hbut], ...
                'min',-1,'max',1,'sliderstep',[1 1],'value',0, ...
                'callback',@(u,e)data_move(X,u));
            X.controls(i).close = uicontrol('parent',hpi,'style','pushbutton', ...
                'units','pixel','position',[W-wplusminus+1+hbut hdec hbut hbut], ...
                'string','X', ...
                'callback',@(u,e)data_remove(X,u));
            display_linecontent(X,i)
        end
        function display_linecontent(X,i)
            % locate the corresponding content and panel
            name = X.x(i).name;
            kspec = find(strcmpi(name,{X.specs.name}));
            if ~isscalar(kspec), error('invalid name'), end
            spec = X.specs(kspec).controls;
            ncontrol = length(spec);
            totallength = sum([spec.length]);
            contentpanel = X.controls(i).contentpanel;
            content = cell(1,ncontrol);
            % sizes
            [W H h hdec hbut hdectext htext] = display_sizes(X);
            pos = get(contentpanel,'position');
            WP = pos(3); 
            wdec = 2;
            wunit = (WP-(ncontrol+1)*wdec)/totallength;
            % display controls
            xpos = wdec+1;
            for k=1:ncontrol
                speck = spec(k);
                if isfield(speck,'label') && ~isempty(speck.label)
                    w = speck.labellength*wunit-1;
                    uicontrol('parent',contentpanel,'style','text', ...
                        'units','pixel','position',[xpos hdectext w htext], ...
                        'string',speck.label);
                    xpos = xpos + w + 1;
                    w = (speck.length-speck.labellength)*wunit;
                else
                    w = speck.length*wunit;
                end
                switch speck.style
                    case {'popupmenu' 'edit' 'pushbutton' 'checkbox'}
                        content{k} = uicontrol('parent',contentpanel, ...
                            'style',speck.style);
                        if isfield(speck,'string'), set(content{k},'string',speck.string); end
                    case 'stepper'
                        content{k} = fn_stepper('parent',contentpanel);
                    otherwise
                        error('style ''%s'' is not handled yet',speck.style)
                end
                set(content{k}, ...
                    'units','pixel','position',[xpos hdec w hbut], ...
                    'callback',@(u,e)data_userset(X,u,k));
                if isfield(speck,'more') && ~isempty(speck.more), set(content{k},speck.more{:}), end
                xpos = xpos + w + wdec;
            end
            X.controls(i).content = content;
            display_contentvalue(X,i)
        end
        function display_contentvalue(X,i)
             % locate the corresponding content and panel
            name = X.x(i).name;
            kspec = find(strcmpi(name,{X.specs.name}));
            if ~isscalar(kspec), error('invalid name'), end
            spec = X.specs(kspec).controls;
            ncontrol = length(spec);
            content = X.controls(i).content;
            for k=1:ncontrol
                speck = spec(k);
                val = X.x(i).value{k};
                switch speck.style
                    case 'popupmenu'
                        klist = find(strcmp(val,speck.string));
                        if ~isscalar(klist), error('value is not in list'), end
                        set(content{k},'value',klist)
                    case 'edit'
                        if ischar(val)
                            set(content{k},'string',val)
                        elseif isnumeric(val)
                            set(content{k},'string','userdata','userdata',val)
                        else
                            error('incorrect value type for edit control')
                        end
                    case 'checkbox'
                        set(content{k},'value',val)
                    case 'pushbutton'
                        set(content{k},'userdata',val);
                    case 'stepper'
                        set(content{k},'value',val)
                    otherwise
                        error('this control style is not handled yet')
                end
            end
        end
        function display_clear(X,ind)
            if nargin<2, ind = 1:length(X.controls); end
            for i = ind
                delete(X.controls(i).panel)
                X.controls(i).panel = [];
            end
        end
        function display_moveline(X,i,j,doswitch)
            [W H h] = display_sizes(X);
            if doswitch
                X.controls([j i]) = X.controls([i j]);
                set(X.controls(i).panel,'position',[1 H-(i+1)*h W h],'userdata',i)
                set(X.controls(j).panel,'position',[1 H-(j+1)*h W h],'userdata',j)
            else
                X.controls(j) = X.controls(i);
                set(X.controls(j).panel,'position',[1 H-(j+1)*h W h],'userdata',j)
            end
        end
    end
    
    % Data
    methods
        function data_userset(X,u,k)
            % function data_userset(X,u,k)
            %---
            % u is the control whose value has been changed by user
            % k is the index of the control in the set of controls 
            i = get(get(get(u,'parent'),'parent'),'userdata'); % index in X.x being modified
            % locate the corresponding content and panel
            name = X.x(i).name;
            kspec = strcmpi(name,{X.specs.name});
            speck = X.specs(kspec).controls(k);
            switch speck.style
                case 'pushbutton'
                    fun = speck.callback;
                    if isempty(fun)
                        disp 'no action for push button'
                        return
                    else
                        vali = X.x_data(i).value;
                        val = feval(fun,vali);
                        if isempty(val), return, end
                    end
                case 'popupmenu'
                    str = speck.string;
                    val = str{get(u,'value')};
                case 'edit'
                    str = get(u,'string');
                    if strcmp(str,'userdata')
                        errordlg('what did you do!?');
                        return
                    else
                        val = str;
                    end
                case 'checkbox'
                    val = logical(get(u,'value'));
                case 'stepper'
                    val = get(u,'value');
                otherwise
                    error('this control style is not handled yet')
            end
            X.x_data(i).value{k} = val;
            % eval fn_supercontrol callback
            X.activechg = X.x(i).active;
            if ~isempty(X.callback), feval(X.callback,X.x), end
        end
        function data_toggleactive(X,u)
            i = get(get(u,'parent'),'userdata');
            X.x_data(i).active = logical(get(u,'value'));
            % eval fn_supercontrol callback
            X.activechg = true;
            if ~isempty(X.callback), feval(X.callback,X.x), end
        end
        function data_addline(X)
            k = get(X.addlinecontrol,'value')-1;
            set(X.addlinecontrol,'value',1)
            if k==0, return, end
            i = X.nx+1;
            X.x_data(i) = struct( ...
                'name',lower(X.specs(k).name), ...
                'active',false, ...
                'value',{{X.specs(k).controls.default}});
            display_line(X,i)
            % eval fn_supercontrol callback
            X.activechg = false;
            if ~isempty(X.callback), feval(X.callback,X.x), end
        end
        function data_move(X,u)
            i = get(get(u,'parent'),'userdata');
            j = i-get(u,'value');
            set(u,'value',0)
            if j<1 || j>X.nx, return, end
            X.x_data([i j]) = X.x_data([j i]);
            display_moveline(X,i,j,true)
            % eval fn_supercontrol callback
            X.activechg = X.x(i).active & X.x(j).active;
            if ~isempty(X.callback), feval(X.callback,X.x), end
        end
        function data_remove(X,u)
            i = get(get(u,'parent'),'userdata');
            X.activechg = X.x(i).active;
            X.x_data(i) = [];
            display_clear(X,i)
            for j=i:X.nx
                display_moveline(X,j+1,j,false)
            end
            % eval fn_supercontrol callback
            if ~isempty(X.callback), feval(X.callback,X.x), end
        end
    end
    
end