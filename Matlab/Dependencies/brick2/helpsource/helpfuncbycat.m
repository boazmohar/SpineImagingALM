%% Functions by Categories
% Brick Toolbox, Version 1.1, 04-Mar-2012
%


%% Array manipulation
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr><td colspan=2><b> indices manipulation</b></tr>
% <tr valign="top"><td width="150"><a href="help_fn_indices.html">         fn_indices        </a></td><td> Convert between global and per-dimension indices
% <tr valign="top"><td width="150"><a href="help_fn_mod.html">             fn_mod            </a></td><td> Return modulus between 1 and n instead of between 0 and n-1
% <tr valign="top"><td width="150"><a href="help_fn_reshapepermute.html">  fn_reshapepermute </a></td><td> Combination of reshape and permute in a single function call
% <tr valign="top"><td width="150"><a href="help_fn_interleave.html">      fn_interleave     </a></td><td> Interleave data
% <tr valign="top"><td width="150"><a href="help_fn_add.html">             fn_add            </a></td><td> Add arrays whose dimensions only partially match
% <tr valign="top"><td width="150"><a href="help_fn_mult.html">            fn_mult           </a></td><td> Multiply arrays whose dimensions only partially match
% <tr valign="top"><td width="150"><a href="help_fn_sizecompare.html">     fn_sizecompare    </a></td><td> Check whether two size vectors are equivalent
% <tr><td colspan=2><b> data operation
% <tr valign="top"><td width="150"><a href="help_fn_bin.html">             fn_bin            </a></td><td> Data binning
% <tr valign="top"><td width="150"><a href="help_fn_round.html">           fn_round          </a></td><td> Round x to the closest multiple of y
% <tr valign="top"><td width="150"><a href="help_fn_coerce.html">          fn_coerce         </a></td><td> Restrict data to a specific range
% <tr valign="top"><td width="150"><a href="help_fn_clip.html">            fn_clip           </a></td><td> Rescale data, restrict the range, color
% <tr valign="top"><td width="150"><a href="help_fn_normalize.html">       fn_normalize      </a></td><td> Normalize data based on averages in some specific dimension
% <tr valign="top"><td width="150"><a href="help_fn_smooth.html">          fn_smooth         </a></td><td> 1D or 2D smoothing using Gaussian convolution
% <tr valign="top"><td width="150"><a href="help_fn_smooth3.html">         fn_smooth3        </a></td><td> 3D smoothing using Gaussian convolution
% <tr><td colspan=2><b> interpolation
% <tr valign="top"><td width="150"><a href="help_fn_decale.html">          fn_decale         </a></td><td> Interpolate a translated vector
% <tr valign="top"><td width="150"><a href="help_fn_translate.html">       fn_translate      </a></td><td> Interpolate a translated image or movie
% <tr valign="top"><td width="150"><a href="help_fn_interprows.html">      fn_interprows     </a></td><td> Realign temporally the frames of a raster scan movie
% <tr valign="top"><td width="150"><a href="help_fn_gettrend.html">        fn_gettrend       </a></td><td> Remove signals modelled by specific low-frequency regressors
% <tr valign="top"><td width="150"><a href="help_fn_register.html">        fn_register       </a></td><td> Coregister frames of a movie
% <tr><td colspan=2><b> min, max, mean,...
% <tr valign="top"><td width="150"><a href="help_fn_max.html">             fn_max            </a></td><td> Global maximum of an array and its coordinates
% <tr valign="top"><td width="150"><a href="help_fn_min.html">             fn_min            </a></td><td> Global minimum of an array and its coordinates
% <tr valign="top"><td width="150"><a href="help_fn_minmax.html">          fn_minmax         </a></td><td> Basic min/max operations (e.g. intersection of 2 ranges)
% <tr valign="top"><td width="150"><a href="help_fn_localmax.html">        fn_localmax       </a></td><td> Find local maxima in a vector
% <tr valign="top"><td width="150"><a href="help_fn_means.html">           fn_means          </a></td><td> Average the successive argumments
% <tr valign="top"><td width="150"><a href="help_fn_meanc.html">           fn_meanc          </a></td><td> Return mean and alpha level confidence interval
% <tr valign="top"><td width="150"><a href="help_fn_meanangle.html">       fn_meanangle      </a></td><td> Average of angles (result in [-pi pi])
% <tr valign="top"><td width="150"><a href="help_fn_meantc.html">          fn_meantc         </a></td><td> Average time course of a movie
% <tr valign="top"><td width="150"><a href="help_fn_trigger.html">         fn_trigger        </a></td><td> Trigger-averaging of movie
% <tr><td colspan=2><b> other
% <tr valign="top"><td width="150"><a href="help_fn_sym.html">             fn_sym            </a></td><td> Convert a symmetric matrix to a vector and vice-versa
% <tr valign="top"><td width="150"><a href="help_fn_timevector.html">      fn_timevector     </a></td><td> Convert set of times to vector of counts and vice-versa
% </table>
% </html>

%% Matlab types
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr><td colspan=2><b> conversions
% <tr valign="top"><td width="150"><a href="help_fn_num2str.html">         fn_num2str        </a></td><td> Convert numeric to char, unless input is already char!
% <tr valign="top"><td width="150"><a href="help_fn_str2double.html">      fn_str2double     </a></td><td> Convert to numeric if not already numeric
% <tr valign="top"><td width="150"><a href="help_fn_str2struct.html">      fn_str2struct     </a></td><td> Evaluate strings that define a structure
% <tr valign="top"><td width="150"><a href="help_fn_struct2str.html">      fn_struct2str     </a></td><td> Convert a structure to strings that define it
% <tr><td colspan=2><b> structures
% <tr valign="top"><td width="150"><a href="help_fn_structdisp.html">      fn_structdisp     </a></td><td> Recursive display of a structure content
% <tr valign="top"><td width="150"><a href="help_fn_structedit.html">      fn_structedit     </a></td><td> Edit a structure; this is a wrapper of fn_control function
% <tr valign="top"><td width="150"><a href="help_fn_structexplore.html">   fn_structexplore  </a></td><td> Navigate using command line inside a large structure
% <tr valign="top"><td width="150"><a href="help_fn_structmerge.html">     fn_structmerge    </a></td><td> Merge two structure
% </table>
% </html>


%% Mathematics
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr><td colspan=2><b> shortcuts
% <tr valign="top"><td width="150"><a href="help_fn_fit.html">             fn_fit            </a></td><td> Shortcut for using Matlab fit function
% <tr><td colspan=2><b> tools
% <tr valign="top"><td width="150"><a href="help_fn_fftfrequencies.html">  fn_fftfrequencies </a></td><td> Frequencies corresponding to the output of Matlab fft function
% </table>
% </html>


%% Programming
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr><td colspan=2><b> shortcuts
% <tr valign="top"><td width="150"><a href="help_fn_switch.html">          fn_switch         </a></td><td> Shortcut for avoiding using if/else and switch
% <tr valign="top"><td width="150"><a href="help_fn_isemptyc.html">        fn_isemptyc       </a></td><td> Which elements of a cell array are empty
% <tr valign="top"><td width="150"><a href="help_fn_disp.html">            fn_disp           </a></td><td> Display multiple arguments
% <tr valign="top"><td width="150"><a href="help_fn_dispandexec.html">     fn_dispandexec    </a></td><td> Display commands in Matlab and executes them
% <tr valign="top"><td width="150"><a href="help_fn_subsref.html">         fn_subsref        </a></td><td> Shortcut for calling Matlab subsref function
% <tr valign="top"><td width="150"><a href="help_fn_ismemberstr.html">     fn_ismemberstr    </a></td><td> Check whether string is part of a set of strings
% <tr valign="top"><td width="150"><a href="help_fn_flags.html">           fn_flags          </a></td><td> Detect flags in the arguments of a function
% <tr valign="top"><td width="150"><a href="help_fn_map.html">             fn_map            </a></td><td> Apply a given function to the elements, columns or rows of an array
% <tr><td colspan=2><b> tools
% <tr valign="top"><td width="150"><a href="help_fn_progress.html">        fn_progress       </a></td><td> Print the state of a calculation
% <tr valign="top"><td width="150"><a href="help_fn_hash.html">            fn_hash           </a></td><td> Unique hash number for an array/cell/structure (Copyright M Kleder)
% <tr><td colspan=2><b> debugging
% <tr valign="top"><td width="150"><a href="help_fn_dbstack.html">         fn_dbstack        </a></td><td> Display current function name, with indent according to stack length
% <tr valign="top"><td width="150"><a href="help_fn_basevars.html">        fn_basevars       </a></td><td> Load base workspace variables in caller workspace and vice-versa
% </table>
% </html>


%% Files
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr><td colspan=2><b> shortcuts
% <tr valign="top"><td width="150"><a href="help_fn_cd.html">              fn_cd             </a></td><td> User definition of shortcut to fast access directories
% <tr valign="top"><td width="150"><a href="help_fn_fileparts.html">       fn_fileparts      </a></td><td> Get specific file parts
% <tr valign="top"><td width="150"><a href="help_fn_ls.html">              fn_ls             </a></td><td> Return folder content
% <tr valign="top"><td width="150"><a href="help_fn_mkdir.html">           fn_mkdir          </a></td><td> Create a directory if does not exist
% <tr valign="top"><td width="150"><a href="help_fn_movefile.html">        fn_movefile       </a></td><td> Rename files in current directory using regular expression
% <tr><td colspan=2><b> user selection
% <tr valign="top"><td width="150"><a href="help_fn_getfile.html">         fn_getfile        </a></td><td> Select file and remember the containing folder of the last selected file
% <tr valign="top"><td width="150"><a href="help_fn_savefile.html">        fn_savefile       </a></td><td> User select file for saving and remember last containing folder
% <tr valign="top"><td width="150"><a href="help_fn_getdir.html">          fn_getdir         </a></td><td> Select directory and remember last containing folder
% <tr><td colspan=2><b> input
% <tr valign="top"><td width="150"><a href="help_fn_readtext.html">        fn_readtext       </a></td><td> Read text file
% <tr valign="top"><td width="150"><a href="help_fn_readasciimatrix.html"> fn_readasciimatrix</a></td><td> Read 2D array from text file
% <tr valign="top"><td width="150"><a href="help_fn_readbin.html">         fn_readbin        </a></td><td> Read binary file containing some header followed by numerical data
% <tr valign="top"><td width="150"><a href="help_fn_readdatlabview.html">  fn_readdatlabview </a></td><td> Read binary matrix (Labview format: nrow, ncol and then data)
% <tr valign="top"><td width="150"><a href="help_fn_readimg.html">         fn_readimg        </a></td><td> Read image
% <tr valign="top"><td width="150"><a href="help_fn_readmovie.html">       fn_readmovie      </a></td><td> Read AVI movie
% <tr valign="top"><td width="150"><a href="help_fn_readmesh.html">        fn_readmesh       </a></td><td> Read surface mesh from file (Anatomist format)
% <tr valign="top"><td width="150"><a href="help_fn_readtexture.html">     fn_readtexture    </a></td><td> Read surface texture from file (Anatomist format)
% <tr valign="top"><td width="150"><a href="help_fn_readxml.html">         fn_readxml        </a></td><td> Read XML file (Copyright Jarek Tuszynski)
% <tr valign="top"><td width="150"><a href="help_fn_readxmllabview.html">  fn_readxmllabview </a></td><td> Read Labview data saved in XML format
% <tr><td colspan=2><b> output
% <tr valign="top"><td width="150"><a href="help_fn_savetext.html">        fn_savetext       </a></td><td> Save text in text file
% <tr valign="top"><td width="150"><a href="help_fn_saveasciimatrix">      fn_saveasciimatrix</a></td><td> Save 2D array in text file
% <tr valign="top"><td width="150"><a href="help_fn_saveimg.html">         fn_saveimg        </a></td><td> Save image or stack of images, options for clipping, color map...
% <tr valign="top"><td width="150"><a href="help_fn_savemovie.html">       fn_savemovie      </a></td><td> Save a movie into an AVI file
% <tr valign="top"><td width="150"><a href="help_fn_savemesh.html">        fn_savemesh       </a></td><td> Save surface as a mesh (Anatomist format)
% <tr valign="top"><td width="150"><a href="help_fn_savetexture.html">     fn_savetexture    </a></td><td> Save surface texture (Anatomist format)
% <tr valign="top"><td width="150"><a href="help_fn_savexml.html">         fn_savexml        </a></td><td> Save a structure in an XML file (Copyright Jarek Tuszynski)
% <tr valign="top"><td width="150"><a href="help_fn_savefig.html">         fn_savefig        </a></td><td> Save figure with possibility to control the output size
% <tr valign="top"><td width="150"><a href="help_linux_savefig.html">      linux_savefig     </a></td><td> Save a figure exactly as it is by calling system operator function
% </table>
% </html>


%% Image operation
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr><td colspan=2><b> basic operations
% <tr valign="top"><td width="150"><a href="help_fn_imvect.html">          fn_imvect         </a></td><td> Convert an image to a vector of pixels inside a mask, and vice-versa
% <tr valign="top"><td width="150"><a href="help_fn_imageop.html">         fn_imageop        </a></td><td> Apply a series of transformations to an image
% <tr valign="top"><td width="150"><a href="help_fn_maskavg.html">         fn_maskavg        </a></td><td> Cluster pixels in an image or movie
% <tr><td colspan=2><b> GUI programs
% <tr valign="top"><td width="150"><a href="help_fn_maskselect.html">      fn_maskselect     </a></td><td> Manual selection of a mask inside an image
% <tr valign="top"><td width="150"><a href="help_fn_subrect.html">         fn_subrect        </a></td><td> Manual selection of a rectangular mask inside an image
% <tr valign="top"><td width="150"><a href="help_fn_alignimage.html">      fn_alignimage     </a></td><td> Manual alignment of 2 images
% <tr valign="top"><td width="150"><a href="help_fn_color2bw.html">        fn_color2bw       </a></td><td> Reduce color to grayscale image while keeping as much information
% </table>
% </html>


%% Data display
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr><td colspan=2><b> shortcuts
% <tr valign="top"><td width="150"><a href="help_fn_drawpoly.html">        fn_drawpoly       </a></td><td> Shortcut for line(poly(:,1),poly(:,2))
% <tr valign="top"><td width="150"><a href="help_fn_isfigurehandle.html">  fn_isfigurehandle </a></td><td> Is handle a plausible figure handle
% <tr><td colspan=2><b> general tools
% <tr valign="top"><td width="150"><a href="help_fn_colorset.html">        fn_colorset       </a></td><td> Set of color, alternative to Matlab default 'ColorOrder'
% <tr valign="top"><td width="150"><a href="help_fn_lines.html">           fn_lines          </a></td><td> Draw a series of vertical and/or horizontal lines
% <tr valign="top"><td width="150"><a href="help_fn_review.html">          fn_review         </a></td><td> Navigate with arrow keys inside a set of data
% <tr valign="top"><td width="150"><a href="help_fn_review_showres.html">  fn_review_showres </a></td><td> function called by fn_review; can be edited by user
% <tr><td colspan=2><b> time courses displays
% <tr valign="top"><td width="150"><a href="help_fn_errorbar.html">        fn_errorbar       </a></td><td> Display nice error bars
% <tr valign="top"><td width="150"><a href="help_fn_drawspline.html">      fn_drawspline     </a></td><td> Fit a curve using splines with movable control points
% <tr valign="top"><td width="150"><a href="help_fn_eegplot.html">         fn_eegplot        </a></td><td> Display multiple time courses spaced from each other
% <tr valign="top"><td width="150"><a href="help_fn_eegdisplay.html">      fn_eegdisplay     </a></td><td> Joint image and time courses display of 2D data
% <tr valign="top"><td width="150"><a href="help_fn_regression.html">      fn_regression     </a></td><td> Display of data points together with linear regression
% <tr valign="top"><td width="150"><a href="help_fn_rasterplot.html">      fn_rasterplot     </a></td><td> Raster plot (punctual events shown as small bars)
% <tr><td colspan=2><b> time courses tools
% <tr valign="top"><td width="150"><a href="help_fn_axis.html">            fn_axis           </a></td><td> Set axis range for best visual aspect
% <tr valign="top"><td width="150"><a href="help_fn_nicegraph.html">       fn_nicegraph      </a></td><td> Improve aspect of graph display
% <tr valign="top"><td width="150"><a href="help_fn_labels.html">          fn_labels         </a></td><td> Improve aspect of graph and add labels and more
% <tr valign="top"><td width="150"><a href="help_fn_plotscale.html">       fn_plotscale      </a></td><td> Two-directional scale bar for graph display
% <tr valign="top"><td width="150"><a href="help_fn_linespecs.html">       fn_linespecs      </a></td><td> Handle abbreviated plot options (e.g. 'r.')
% <tr><td colspan=2><b> 2D displays
% <tr valign="top"><td width="150"><a href="help_fn_displayarrows.html">   fn_displayarrows  </a></td><td> Display an image and a velocity field on top of it
% <tr valign="top"><td width="150"><a href="help_fn_tensordisplay.html">   fn_tensordisplay  </a></td><td> Display a field of 2x2 symmetric matrices using ellipses
% <tr valign="top"><td width="150"><a href="help_fn_framedisplay.html">    fn_framedisplay   </a></td><td> Sequential display of frames from a movie
% <tr><td colspan=2><b> 2D tools
% <tr valign="top"><td width="150"><a href="help_fn_axispixel.html">       fn_axispixel      </a></td><td> Set axis size for an optimal display of images
% <tr valign="top"><td width="150"><a href="help_fn_imdistline.html">      fn_imdistline     </a></td><td> Show the distance between two points (Copyright The MathWorks)
% <tr valign="top"><td width="150"><a href="help_fn_scale.html">           fn_scale          </a></td><td> Scale bar for image display
% <tr valign="top"><td width="150"><a href="help_fn_showcolormap.html">    fn_showcolormap   </a></td><td> Display a color map
% <tr><td colspan=2><b> movie displays
% <tr valign="top"><td width="150"><a href="help_fn_playmovie.html">       fn_playmovie      </a></td><td> Simple showing of a movie
% <tr valign="top"><td width="150"><a href="help_fn_movie.html">           fn_movie          </a></td><td> Show a movie, large number of options
% <tr><td colspan=2><b> mesh computations and displays
% <tr valign="top"><td width="150"><a href="help_fn_meshclosestpoint.html">fn_meshclosestpoint</a></td><td> Closest vertex on a mesh to a given point
% <tr valign="top"><td width="150"><a href="help_fn_meshinv.html">         fn_meshinv        </a></td><td> Invert the orientation of triangular faces
% <tr valign="top"><td width="150"><a href="help_fn_meshnormals.html">     fn_meshnormals    </a></td><td> Compute the normals to faces of a mesh
% <tr valign="top"><td width="150"><a href="help_fn_meshplot.html">        fn_meshplot       </a></td><td> Display a mesh
% <tr valign="top"><td width="150"><a href="help_fn_meshselectpoint.html"> fn_meshselectpoint</a></td><td> Display a mesh and let user select a point with mouse
% <tr valign="top"><td width="150"><a href="help_fn_cubemesh.html">        fn_cubemesh       </a></td><td> Render the "faces" of a 3D data (creates a mesh and texture)
% <tr valign="top"><td width="150"><a href="help_fn_cubeview.html">        fn_cubeview       </a></td><td> Render the "faces" of a 3D data (creates an image)
% <tr><td colspan=2><b> elaborate programs
% <tr valign="top"><td width="150"><a href="help_fn_imvalue.html">         fn_imvalue        </a></td><td> Automatic link graphs and images for point selection and zooming
% <tr valign="top"><td width="150"><a href="help_fn_4Dview.html">          fn_4Dview         </a></td><td> Navigation inside 3D, 4D or 5D imaging data
% </table>
% </html>


%% GUI programming
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr><td colspan=2><b> shortcuts
% <tr valign="top"><td width="150"><a href="help_fn_evalcallback.html">    fn_evalcallback   </a></td><td> Evaluate a callback, i.e. a char array, function handle or cell array
% <tr valign="top"><td width="150"><a href="help_fn_multcallback.html">    fn_multcallback   </a></td><td> Apply several callback functions to an object
% <tr valign="top"><td width="150"><a href="help_fn_get.html">             fn_get            </a></td><td> Get mutiple properties of multiple objects at once
% <tr valign="top"><td width="150"><a href="help_fn_set.html">             fn_set            </a></td><td> Set simultaneously multiple properties of multiple graphic objects
% <tr><td colspan=2><b> tools
% <tr valign="top"><td width="150"><a href="help_fn_pixelpos.html">            fn_pixelpos           </a></td><td> Position of an object in pixel units
% <tr valign="top"><td width="150"><a href="help_fn_pixelsize.html">           fn_pixelsize          </a></td><td> Size of an object in pixel units
% <tr valign="top"><td width="150"><a href="help_fn_coordinates.html">         fn_coordinates        </a></td><td> Conversion of screen/figure/axes, normalized/pixel coordinates
% <tr valign="top"><td width="150"><a href="help_fn_controlpositions.html">    fn_controlpositions   </a></td><td> Position uicontrol objects relative to an axes
% <tr valign="top"><td width="150"><a href="help_fn_setfigsize.html">          fn_setfigsize         </a></td><td> Change the size of a figure and check that it fits the screen
% <tr valign="top"><td width="150"><a href="help_fn_setpropertyandmark.html">  fn_setpropertyandmark </a></td><td> Change both an object property and a control value
% <tr><td colspan=2><b> mouse actions
% <tr valign="top"><td width="150"><a href="help_fn_buttonmotion.html">    fn_buttonmotion   </a></td><td> Execute a task while mouse pointer is moved around
% <tr valign="top"><td width="150"><a href="help_fn_moveobject.html">      fn_moveobject     </a></td><td> Move a graphic object with mouse
% <tr valign="top"><td width="150"><a href="help_fn_getline.html">         fn_getline        </a></td><td> Manual selection of a polygon (Copyright The MathWorks)
% <tr valign="top"><td width="150"><a href="help_fn_mouse.html">           fn_mouse          </a></td><td> Manual selection of a variety of shapes
% <tr><td colspan=2><b> elaborate tools
% <tr valign="top"><td width="150"><a href="help_fn_framedesign.html">     fn_framedesign    </a></td><td> Utility to let user reposition graphic objects inside a figure
% <tr valign="top"><td width="150"><a href="help_interface.html">          interface         </a></td><td> Parent class to create cool graphic interfaces
% <tr><td colspan=2><b> pre-defined arrangements of controls
% <tr valign="top"><td width="150"><a href="help_fn_okbutton.html">        fn_okbutton       </a></td><td> Small 'ok' button waits to be pressed
% <tr valign="top"><td width="150"><a href="help_fn_menu.html">            fn_menu           </a></td><td> Utility to create a basic GUI made of a line of buttons
% <tr><td colspan=2><b> special controls
% <tr valign="top"><td width="150"><a href="help_fn_multcheck.html">       fn_multcheck      </a></td><td> Special control made of multiple check boxes
% <tr valign="top"><td width="150"><a href="help_fn_buttongroup.html">     fn_buttongroup    </a></td><td> Set of radio buttons or toggle buttons
% <tr valign="top"><td width="150"><a href="help_fn_slider.html">          fn_slider         </a></td><td> Special control that improves the functionality of Matlab slider
% <tr valign="top"><td width="150"><a href="help_fn_sliderenhance.html">   fn_sliderenhance  </a></td><td> Allow a slider uicontrol to evaluate its callback during scrolling
% <tr valign="top"><td width="150"><a href="help_fn_stepper.html">         fn_stepper        </a></td><td> Special numeric control that includes increment/decrement buttons
% <tr valign="top"><td width="150"><a href="help_fn_sensor.html">          fn_sensor         </a></td><td> Special control whose value is changed by moving the mouse
% <tr><td colspan=2><b> elaborate controls
% <tr valign="top"><td width="150"><a href="help_fn_control.html">         fn_control        </a></td><td> Arrangement of controls that reflect the state of a set of parameters
% <tr valign="top"><td width="150"><a href="help_fn_supercontrol.html">    fn_supercontrol   </a></td><td> Super arrangement of controls
% <tr valign="top"><td width="150"><a href="help_fn_reallydlg.html">       fn_reallydlg      </a></td><td> Ask for confirmtaion
% <tr valign="top"><td width="150"><a href="help_fn_dialog_questandmem.html"> fn_dialog_questandmem </a></td><td> Confirmation question with an option for not asking again
% <tr valign="top"><td width="150"><a href="help_fn_input.html">           fn_input          </a></td><td> Prompt user for a single value. Function based on fn_structedit
% </table>
% </html>


%% Miscellaneous
% <html>
% <table cellspacing="0" width="100%" border="0" cellpadding="2" style="margin-top:5px;margin-bottom:30px;">
% <tr valign="top"><td width="150"><a href="help_alias.html">              alias             </a></td><td> Create command shortcuts
% <tr valign="top"><td width="150"><a href="help_fn_email.html">           fn_email          </a></td><td> Send e-mails from Matlab! Can attach figures, M-files and more
% <tr valign="top"><td width="150"><a href="help_fn_figmenu.html">         fn_figmenu        </a></td><td> An automatic custom menu for figures: save figure, distance tool, ...
% <tr valign="top"><td width="150"><a href="help_pointer.html">            pointer           </a></td><td> Implement a pointer to any Matlab object
% </table>
% </html>



