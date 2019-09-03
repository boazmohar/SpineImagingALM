%% fn_4Dview

%% Syntax
%  fn_4Dview(option1,optionArgs1,...)
%  varargout = fn_4Dview(action,actionArgs)

%% Description
%  displays spatial or temporal slices of multi-dimensional data
%  by default: data is considered to be 3D, 4D or 5D (dimensions ordered as
%  x-y-z-t-m) and is displayed as 3 spatial sections at a given time point
%  use option flags for other display types
% 
%  type 'fn_4Dview demo' to see a demo
% 
%  Input:
%  - option..        string - specify an option (see below)
%  - optionArgs..    argument of option
%  - action..        specific usage of fn_4Dview, does not create a new
%                    interface (see below)
% 
%  Options:
%  DATA
%  - 'data',data     array, by default, dimensions should be x - y - z - t - m 
%                    (m is used for multiple data or channels)   
%                    [the 'data' flag can be omitted]
%                    [default: data=0 or zeros(appropriate size)]
%  - 'options',{options}     some display types allow additional options like 
%                    'xxxplot' -> same option possibilities as the Matlab plot function, 
%                    'quiver'...     
%                    [the 'options' flag can be omitted] [default: {}]
%  TYPE OF DISPLAY
%  - '3d'            data is three-dimensional (dimensions are x - y - z - t - m)
%                    and will be spatially displayed as three cross
%                    sections, for a given time point
%                    [default]
%  - '3dplot'        data is three-dimensional (dimensions are x - y - z - t - m)
%                    and will be temporally displayed as a time course, for
%                    a given space point
%  - '2d'            data is two-dimensional (dimensions are x - y - t - m)
%                    and will be spatially displayed as a single image
%  - '2dplot'        data is two-dimensional (dimensions are x - y - t - m)
%                    and will be displayed temporally
%  - '0d',tidx       [deprecated]
%  - 'quiver'        data is a 2D vector field, plus possibly 2D image
%                    (dimensions are x - y - t - [2 vector components +
%                    image component])
%                    this is only applicable for spatial display
%  - 'mesh',mesh     data is one-dimensional (dimensions are i - t - m),
%                    where i refers to the ith vertex of a 3D mesh specified
%                    as a cell array or struct following the 'mesh' flag 
%                    {[3 x nv vertices],[3 x nt triangles]},
%                    and will be spatially displayed as a surface
%  - 'meshplot',mesh data is one-dimensional (dimensions are i - t - m) as
%                    for 'mesh', and will be displayed temporally
%  - 'ind',indices   data is one-dimensional (dimensions are i - t - m),  
%                    indices is a 2D or 3D image which entails indices of
%                    the first data component (see example)
%                    this is only applicable for temporal display
%  - 'timeslider',tidx   
%                    creates a slider control to move in time array tidx
%  - 'ext',{@updatefcn,par1,par2,...}
%    'ext',{'command','infoname'}
%                    fn_4Dview does not display anything, but it links the
%                    execution of function updatefcn to other objects
%                    handled by fn_4Dview. The function prototype should be:
%                        updatefcn = function(info,par1,par2,...)
%                    with info being the structure described in fn_4Dview code.
%                    Alternatively, one can use a string command; each time
%                    before the command is executed, info is stored in the
%                    base workspace with the name 'infoname', so that the
%                    command can access to the informations it encloses.
%                    Attention, the information concerning these updates
%                    needs to be stored in a graphical object; if none is
%                    specified (using 'in'), one is automatically created,
%                    but one should delete it later on using 'fn_4Dview
%                    unregister' (see below)
%  BEHAVIOR
%  - 'active'        allows callbacks (e.g. selecting point with mouse) [default]
%  - 'passive'       disallows callbacks
%  - 'key',k         use this to have independent links between windows
%                    [default: all windows are linked using same key k=1]
%  DATA TRANSFORM
%  - 'mat',M         defines a spatial linear transformation between
%                    data and world coordinates [default = eye(4)] (pixel
%                    indices start at 1)
%                    can be a 4x4 or 3x4 matrix (rotation+translation) or
%                    3x3 matrix (rotation only) or 3x1 vector (translation
%                    only) or 1x3,1x2,1x1 vector (special: scaling, first
%                    pixel at [0 0 0])
%  - 'tidx',[t0 dt]  time of first frame and interval between successive
%                    frames [default: t0-=1, dt=1]
%  - 'tidx',tidx     alternatively, one can define every sampling points 
%                    time coordinates (tidx should not have 2 elements);
%                    it is necessary that they are equally spaced then
%  - 'heeginv',H     provides a matrix to multiply data with before using it
%                    this is only applicable for spatial mesh display
%  - 'applyfun',fun  function handler or {function handler, additional arguments, ...}
%                    data will be transformed according to fun before every
%                    display - NOT IMPLEMENTED YET
%  DISPLAY
%  - 'in',handle     forces display in figure, axes or uicontrol specified by handle
%                    [the 'in' flag can be ommited] 
%                    [default: in active figure or axes]
%  - 'channel',f       if data has a 'multiple data' component (dimension m is
%                    non-singleton), specifies which one to use for spatial
%                    display [default: f=1]
%  - 'clip',[m M]    specify a clipping for image display
%  - 'xyzperm',[a b c] defines a permutation of the dimensions before
%                    visualizing the 3 cross-sections of the data
%  - 'zfact',r       ratio between z and xy resolutions for 3D display
%  - 'labels',{xlabel,ylabel,zlabel} axis labels (ylabel and zlabel can be omitted)
% 
%  Actions:
%  - 'demo'                  run the demo
%  - 'changexyz',key,xyz     update space coordinates for objects linked by the key       
%  - 'changet',key,t         update time coordinate for objects linked by the key
%  - 'unregister'[,@updatefcn][,key] 
%                            unregister external functions which were
%                            previously registered (type 'ext') by deleting
%                            the associated objects;
%                            this unregistration can be filtered by which
%                            function handles and/or which linking key

%% Source
% Thomas Deneux
%
% Copyright 2004-2012
%
% last modification: September 20th, 2007
%
