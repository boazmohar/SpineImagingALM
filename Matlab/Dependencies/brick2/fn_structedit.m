function s = fn_structedit(s,varargin)
% function s = fn_structedit(s[,spec][,hp])
% function s = fn_structedit(fieldname1,value1,fieldname2,value2,...)
%---
% allow user to interactively modify a structure, and returns the modified
% structure; returns an empty array if the figure is closed
% this function is a wrapping of fn_control
% 
% See all fn_control, fn_input, fn_reallydlg

% Thomas Deneux
% Copyright 2007-2012

if ischar(s)
    s = struct(s,varargin{:});
    varargin = {};
end
X = fn_control(s,varargin{:},'ok');
addlistener(X,'OK',@(u,e)update)
s = [];
% wait that the figure will be destroyed
waitfor(X.hp)

% nested function: update of s
function update
    s = X.s;
end

end
        
        

