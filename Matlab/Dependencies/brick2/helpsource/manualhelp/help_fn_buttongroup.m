%% fn_buttongroup
% Create an array of radion buttons or toggle buttons.
%
%% Syntax
%  G = fn_buttongroup(style,str,callback,'prop1',value1,...)
%
%% Input:
%
%  - style           'radio' or 'toggle'
%  - str             list of string values
%  - callback        function with prototype @(x)fun(x), where x is the
%                    selected string value
%  - propn/valuen    additional properties to be set (possibilities are:
%                    'parent', 'units', 'position', 'value')
%
%% Example

fn_buttongroup('radio', ...
    {'science' 'France' 'God' 'you' 'fn_buttongroup'}, ...
    @(x)disp(['I love ' x]), ...
    'parent',figure('pos',[400 200 140 180]), ...
    'value',4);

%% Source
% Thomas Deneux
%
% Copyright 2010-2012   
