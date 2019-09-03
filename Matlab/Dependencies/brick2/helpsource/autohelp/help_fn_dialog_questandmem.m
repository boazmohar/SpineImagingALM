%% fn_dialog_questandmem

%% Syntax
%  answer = fn_dialog_questandmem(question[,identifier])

%% Description
%  This imitates Windows style confirmation window, with a choice between
%  'Yes' and 'Cancel', and additionally a possibility to mark a box 'Do not
%  ask again'.
% 
%  Input:
%  - question        string - the question
%  - identifier      string - a unique identifier that serves to store the
%                    choice not to ask again the question [default: the
%                    question is used as an identifier]
% 
%  Output:
%  - answer          logical - true for 'Yes', false for 'Cancel'

%% Source
% Thomas Deneux
%
% Copyright 2010-2012
%
