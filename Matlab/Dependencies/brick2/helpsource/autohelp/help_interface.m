%% interface

%% Description
%  Interface class provide utilities for designing a graphic interface,
%  such as allowing the user to resize the graphical elements, loading
%  and auto-saving program options, etc..
% 
%  Notes:
%  - to make a new interface, define a new class having interface as a
%    parent
%  - constructor: 
%      . in the new object constructor, first call the interface
%        constructor (X = X@interface(hf,figtitle,defaultoptions)
%      . then define the graphic objects that user can resize in the
%        'grob' property
%      . then call interface_end(X) to auto-position these objects
%  - methods:
%      . if the child class defines menus additional to the one of
%        interface, it should do it in a init_menus(X), which starts by
%        calling init_menus@interface(X); menu handles can be stored in
%        the structure X.menus
%      . interface overwrites the default set method in order to easily 
%        provide the user with a description of the value to enter; for
%        such to happen, the child class should have a method x =
%        setinfo(X) that returns a stucture with field names and
%        description values (which can be a string or a cell with
%        possible values)
% 
%  An example is provided in [Brick toolbox dir]/examples/example_interface

%% Source
% Thomas Deneux
%
% Copyright 2007-2012
%
