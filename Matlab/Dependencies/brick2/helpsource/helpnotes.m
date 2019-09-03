%% Useful notes on the Brick Toolbox
%
% We hope the 'Brick' Toolbox will help you, its documentation is not
% complete, so you won't have too much to read! Just pay attention to the
% following points:

%% Higlights
% The <fn_demo.html main presentation> has already highlithed functions
% |<help_fn_coerce.html fn_coerce>|, 
% |<help_fn_add.html fn_add>|,
% |<help_fn_mult.html fn_mult>|,
% |<help_fn_normalize.html fn_normalize>|,
% |<help_fn_bin.html fn_bin>|,
% |<help_fn_flags.html fn_flags>|,
% |<help_fn_switch.html fn_switch>|, 
% |<help_fn_errorbar.html fn_errorbar>|,
% |<help_fn_buttonmotion.html fn_buttonmotion>|,
% |<help_fn_control.html fn_control>|, 
% |<help_fn_imvalue.html fn_imvalue>|, 
% |<help_fn_figmenu.html fn_figmenu>| and
% |<help_fn_email.html fn_email>|. 
%
% Here are a few additional functions that are extremely useful: 
%
% *|<help_fn_getfile.html fn_getfile>|*:
% <matlab:fn_dispandexec('fn_getfile') try it>, choose file in a different
% directory than the one that first pops out, then try again to
% understand the difference with Matlab |uigetfile|
% 
% *|<help_fn_cd.html fn_cd>|*:
% define shorcuts to your favorite folders by typing 
% |<matlab:fn_dispandexec('fn_cd(''edit'')') fn_cd edit>|.
%
% *|<help_fn_progress.html fn_progress>|*:
% try 
% <matlab:fn_dispandexec({'fn_progress(''executed'',100)','for$k=1:100,$pause(.03),$fn_progress(k),$end'})
% this>
%
% *|<help_fn_movie.html fn_movie>|*:
% try the
% <matlab:fn_dispandexec('fn_movie(''demo'')') demo>, press 'speed' and play
% with the controls
%
% *|<help_fn_4Dview.html fn_4Dview>|*: 
% try
% <matlab:fn_dispandexec({'load(fullfile(fileparts(which(''fn_movie'')),''data/fn_movie_demo''))','fn_4Dview(Y,''in'',subplot(121),''2d''),colormap(''gray'')','fn_4Dview(Y,''in'',subplot(122),''2dplot'')'})
% this>
%
% *|<help_interface.html interface>|*:
% This is a parent class that will help creating GUIs, especially regarding
% the positioning of graphic objects (see also function
% |<help_fn_framedesign.html fn_framedesign>|). To see how it works, check the example
% <matlab:edit(fullfile(fileparts(which('fn_movie')),'examples/interface_example'))
% interface_example.m>
% Note that probably the most developed side of the Brick Toolbox is its
% utilies for creating GUIs. 
%
% *|<help_alias.html alias>|*:
% Create your own shortcuts in Matlab. Note that, in this vein, we also
% highly recommend to use the *shortcuts toolbar* in Matlab. Learn about it  
% <"http://blogs.mathworks.com/desktop/2007/03/29/shortcuts-for-commonly-used-code/"
% here>.
%
% Of course, you wil find all functions organized by categories in the
% <helpfuncbycat.html Function Reference>.

%% Installation
% Add only the main Brick folder to the Matlab path (not the
% subdirectories).
%
% A new menu should appear in the Help window (type 'doc'). If it is not the
% case, make sure that in the Matlab Preferences, Help tab, 'All products'
% is selected.

%% Syntax of the functions
% The Brick toolbox does not follow (at least for the moment) the Matlab
% convention of the 'H-line', i.e. that the first line in the help of a
% function should be a short description of what the function is doing.
% Instead, it displays the different syntaxes. Best is to show examples:

%%
%  function b = fn_interprows(a,subrate[,method[,extrap]])

%%
% Brackes |[]| indicate optional arguments. Thus here, |method| and
% |extrap| are optionals; note that it is 
% necessary to provide |method| in order to be able to provide |extrap|

%%
%  function M = fn_savemovie(a[,fname][,clip][,fps][,zoom][,map][,hf])

%%
% Here on the contrary, all the optional can be provided alone (and in most
% functions, they can also be provided in all different orders). The
% function figures out by itself for example that
% |fn_savemovie(data,2,15,'mymovie.avi')| means that 2 is a zoom factor, 15
% a number of frames per second and 'mymovie.avi' a file name.

%% Backward compatibility
% Code that uses this version of the Brick toolbox should be able to run on
% the future versions as well, since the syntax of functions is not supposed to be
% changed (or only enhanced). If nevertheless such changes occur, they will
% be signaled.




