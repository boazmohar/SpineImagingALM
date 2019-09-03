%% fn_buttonmotion
% Executes a task while the mouse is moving, and until
% the mouse button is released. 
%
%% Syntax
%  [varargout =] fn_buttonmotion(fun[,hf])
%  fn_buttonmotion demo


%% Example
% The example below displays the pointer
% position in the figure whenever the mouse is pressed and moved. 
% (<matlab:fn_buttonmotion('demo') execute>).

figure(1), clf
ht = uicontrol('style','text');
fun = ['p=get(1,''currentpoint''); p=p(1,1:2); ' ...
    'set(ht,''pos'',[p 60 20],''string'',num2str(p))'];
set(1,'buttondownfcn',@(u,evnt)fn_buttonmotion(fun))

%% See also
% <"help_fn_moveobject.html" fn_moveobject>

%% Source
% Thomas Deneux
% Copyright 2007-2012
