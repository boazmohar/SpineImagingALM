classdef interface_example < interface
    % This is an example of how to build a GUI using class syntax and
    % using the functionalities of the parent class 'interface'.
    % It lets user define some text and writes it in big in an axes. It
    % also has some menus where one can choose the color.
    
    properties
        % note that some properties, such as the default options or the
        % handles of graphic objects, are stored in properties of the
        % parent class
        x = ';-)';      % the string (already initialized here)
        col             % the color (not initialized yet)
        menuitems       % the menu items
    end
    
    % Initializations
    methods
        function E = interface_example(value)
            % class constructor
            
            % default options
            defaultoptions = struct('color','r');

            % figure
            hf = figure;
            
            % initial call to the parent constructor
            E = E@interface(hf,'My Interface',defaultoptions);
            
            % init graphic objects - no need to define their positions!!!
            E.grob.hu = uicontrol('parent',hf, ...
                'style','edit','string',E.x, ...
                'callback',@(u,e)set(E,'x',get(u,'string')));
            E.grob.ha = axes('parent',hf);
            
            % call parent method 'interface_end' to handle the positions of
            % objects and create the menus
            interface_end(E)
            
            % it is a nice practice to delete the object once the figure is
            % closed
            set(hf,'deletefcn',@(u,e)delete(E))
            
            % set data value
            % note that the default options have been redefined to their
            % last saved value  
            % note also that this will cause an automatic update of the
            % display
            if nargin>1
                E.x = value;
            end
            E.col = E.options.color;
            
        end
        function init_menus(E)
            % initialization of menus: this function will be called
            % automatically during the call to interface_end(E); indeed,
            % interface creates its own menu with several options, and here
            % we can create additional menus
            
            % first put the default interface menu if we want it
            init_menus@interface(E)
            
            % delete the menu if it already exists (this happens when user
            % selects the 'reinit_menus' option)
            if isfield(E.menus,'color') && ishandle(E.menus.color), delete(E.menus.color), end
            
            % a structure where to store the handles of menu options that
            % need to be stored
            E.menuitems = struct;
            
            % the menu
            m = uimenu('parent',E.hf,'label','color');
            E.menus.color = m;
            
            % it childs
            m1 = uimenu(m,'label','choose color');
            uimenu(m,'label','make current color default','callback',@(u,e)makedefaultcol(E))
            uimenu(m,'label','return to default color','callback',@(u,e)backdefaultcol(E))
            
            % sub-childs: choice of color
            E.menuitems.k = uimenu(m1,'label','black'  ,'callback',@(u,e)set(E,'col','k'));
            E.menuitems.b = uimenu(m1,'label','blue'   ,'callback',@(u,e)set(E,'col','b'));
            E.menuitems.r = uimenu(m1,'label','red'    ,'callback',@(u,e)set(E,'col','r'));
            E.menuitems.g = uimenu(m1,'label','green'  ,'callback',@(u,e)set(E,'col','g'));
            E.menuitems.y = uimenu(m1,'label','yellow' ,'callback',@(u,e)set(E,'col','y'));
            if ~isempty(E.col) % during initialization, E.col is not defined yet
                set(E.menuitems.(E.col),'checked','on')
            end
            
            % save the handles
            E.menus.color = m;
        end
    end
    
    % Actions
    methods
        function set.x(E,value)
            % some checks
            if ~ischar(value), error argument, end
            if strcmp(E.x,value), return, end
            % change property value
            E.x = value;
            % update display
            displaytext(E)
        end
        function set.col(E,c)
            % some checks
            if ~ismember(c,'kbrgy'), error argument, end
            oldcol = E.col;
            if oldcol==c, return, end
            % change property value
            E.col = c;
            % update menu marks
            if ~isempty(oldcol), set(E.menuitems.(oldcol),'checked','off'), end %#ok<MCSUP>
            set(E.menuitems.(c),'checked','on')       %#ok<MCSUP>
            % update display
            displaytext(E)
        end
        function makedefaultcol(E)
            % change the options and save it 
            E.options.color = E.col;
            saveoptions(E)
        end
        function backdefaultcol(E)
            E.col = E.options.color;
        end
        function displaytext(E)
            cla(E.grob.ha)
            text(.1,.2,E.x,'color',E.col,'fontsize',50,'fontweight','bold')
        end
    end
    
end
