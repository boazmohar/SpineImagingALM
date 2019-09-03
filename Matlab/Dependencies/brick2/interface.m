classdef interface < hgsetget
    % Interface class provide utilities for designing a graphic interface,
    % such as allowing the user to resize the graphical elements, loading
    % and auto-saving program options, etc..
    %
    % Notes:
    % - to make a new interface, define a new class having interface as a
    %   parent
    % - constructor: 
    %     . in the new object constructor, first call the interface
    %       constructor (X = X@interface(hf,figtitle,defaultoptions)
    %     . then define the graphic objects that user can resize in the
    %       'grob' property
    %     . then call interface_end(X) to auto-position these objects
    % - methods:
    %     . if the child class defines menus additional to the one of
    %       interface, it should do it in a init_menus(X), which starts by
    %       calling init_menus@interface(X); menu handles can be stored in
    %       the structure X.menus
    %     . interface overwrites the default set method in order to easily 
    %       provide the user with a description of the value to enter; for
    %       such to happen, the child class should have a method x =
    %       setinfo(X) that returns a stucture with field names and
    %       description values (which can be a string or a cell with
    %       possible values)
    %
    % An example is provided in [Brick toolbox dir]/examples/example_interface
    
    % Thomas Deneux
    % Copyright 2007-2012
    
    properties (Access='private')
        settings
        figtitle
    end
    properties (SetAccess='protected')
        hf
        grob
        options = struct('positions',struct('name','base','value',[]));
        menus = struct('interface',[]);
    end
    
    methods
        function I = interface(hf,figtitle,defaultoptions)
            if nargin==0
                hf = 827;
                figtitle = 'INTERFACE';
            end
            
            % options
            I.settings = [which(class(I)) 'at'];
            interfaceoptions = I.options; % minimal default options
            loadoptions(I)                % load saved options
            if nargin>=3 && ~isempty(defaultoptions)
                % some new options might have been created (either in
                % 'interface' or in the child class), other might have been
                % removed -> update the saved options
                defaultoptions = fn_structmerge(interfaceoptions,defaultoptions); % default set
                I.options = fn_structmerge(defaultoptions,I.options,'skip');      % update
                saveoptions(I)
            end
            
            % figure
            I.hf = hf;
            I.figtitle = figtitle;
            figure(hf)
            set(hf,'numbertitle','off','name',figtitle)
            set(hf,'resize','off')
            if ~isempty(I.options.positions(1).value), set(hf,'pos',I.options.positions(1).value.hf), end
            clf(hf), set(hf,'menubar','none')
            I.grob.hf = hf;
            
            % end of initialization
            if strcmp(class(I),'interface'), interface_end(I), end
        end
        
        function interface_end(I)
            % end of initialization: call after additional object
            % initializations
            % these initializations should in particular define grob
            
            % interface menu
            init_menus(I)
            
            % position frames
            chgframepositions(I,'set')
            
            % delete object when closing the window
            set(I.hf,'closeRequestFcn',@deleteandclose)
            function deleteandclose(hf,e) %#ok<INUSD>
                delete(I)
                delete(hf)
            end 
            
            % make object available in base workspace
            assignin('base',inputname(1),I)
        end
    end
    
    methods
        function init_menus(I)
            m = I.menus.interface; 
            if ishandle(m), delete(m), end
            m = uimenu('parent',I.hf,'label',I.figtitle);
            I.menus.interface = m;
            
            % positioning
            uimenu(m,'label','Resize frames', ...
                'callback',@(u,evnt)chgframepositions(I,'reset'));
            pos = I.options.positions; npos = length(pos);
            m1 = uimenu(m,'label','Preset positions');
            for k=1:npos
                f = pos(k).name;
                uimenu(m1,'label',f, ...
                    'callback',@(u,evnt)chgframepositions(I,'load',f));
            end
            uimenu(m1,'label','Create new...','separator','on', ...
                'callback',@(u,evnt)chgframepositions(I,'new'));
            m2 = uimenu(m1,'label','Delete');
            for k=1:npos
                f = pos(k).name;
                uimenu(m2,'label',f, ...
                    'callback',@(u,evnt)chgframepositions(I,'delete',f));
            end
            
            % more
            uimenu(m,'label','Edit code','separator','on', ...
                'callback',@(u,evnt)edit(which(class(I))));
            uimenu(m,'label','Reinit menus', ...
                'callback',@(u,evnt)init_menus(I));
            varname = inputname(1);
            uimenu(m,'label','Object in base workspace', ...
                'callback',@(u,evnt)assignin('base',varname,I));
            uimenu(m,'label','Save PNG','separator','on', ...
                'callback',@(u,evnt)fn_savefig(I.hf, ...
                [fn_cd('capture') '/' get(I.hf,'name') '_' datestr(now,'YYmmDDHHMMSS') '.png']));
        end
        function saveoptions(I)
            x = I.options;  %#ok<NASGU>
            save(I.settings,'-STRUCT','x')
        end
        function loadoptions(I)
            if ~exist(I.settings,'file')
                saveoptions(I)
                return
            end
            I.options = load(I.settings);
            % changes from previous version
            if isfield(I.options,'pos')
                pos = I.options.pos;
                I.options = rmfield(I.options,'pos');
                I.options.positions = struct('name','base','value',pos);
                %saveoptions(I)
            end
        end
        function chgframepositions(I,flag,varargin)
            pos = I.options.positions;
            npos = length(pos);
            curpos = pos(1).value;
            doreinitmenus = false;
            switch flag
                case 'reset'
                    pos(1).value = fn_framedesign(I.grob,curpos,true);
                case 'set'
                    pos(1).value = fn_framedesign(I.grob,curpos,false);
                case 'load'
                    f = varargin{1};
                    idx = find(strcmp(f,{pos.name}));
                    pos = pos([idx setdiff(1:npos,idx)]);
                    pos(1).value = fn_framedesign(I.grob,pos(1).value,false);
                    doreinitmenus = true;
                case 'new'
                    name = inputdlg('Name of new position configuration','Enter name',1);
                    if isempty(name) || fn_ismemberstr(name,{pos.name})
                        errordlg('Invalid name (empty or already exists)')
                        return
                    end
                    pos = [struct('name',name,'value',curpos) pos];
                    doreinitmenus = true;
                case 'delete'
                    f = varargin{1};
                    idx = find(strcmp(f,{pos.name}));
                    if idx==1
                        errordlg('Cannot delete current configuration')
                        return
                    end
                    pos(idx) = []; 
                    doreinitmenus = true;
            end
            I.options.positions = pos;
            saveoptions(I)
            if doreinitmenus, init_menus(I), end
        end
        function set(I,f,x)
            if nargin<3
                desc = setinfo(I);
                if nargin<2
                    M = metaclass(I);
                    M = [M.Properties{:}];
                    for k=1:length(M)
                        if ~strcmp(M(k).SetAccess,'public'), continue, end
                        f = M(k).Name;
                        if isfield(desc,f), str = makestr(desc.(f)); else str=[]; end
                        if isempty(str)
                            fprintf('\t%s\n',f)
                        else
                            fprintf('\t%s: %s\n',f,str)
                        end
                    end
                else
                    if isfield(desc,f), disp(makestr(desc.f)), end
                end
            else
                I.(f) = x;
            end
        end
        function x = setinfo(I) %#ok<MANU>
            x = struct;
        end
        
    end
    
end

function desc = makestr(desc)

if isempty(desc)
    desc = '';
elseif iscell(desc)
    desc = [ '[' sprintf(' %s |',desc{:})];
    desc(end) = ']';
end

end







