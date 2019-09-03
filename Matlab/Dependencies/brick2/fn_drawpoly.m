function hl = fn_drawpoly(poly,varargin)
% function hl = fn_drawpoly(poly,varargin)
%---
% shortcut for line(poly(:,1),poly(:,2),varargin{:})

% Thomas Deneux
% Copyright 2006-2012

if nargin==0, help fn_drawpoly, end

if size(poly,1)~=2, error('poly should have two rows'), end
opt = fn_linespecs(varargin{:});
hl = line(poly(1,:),poly(2,:),opt{:});

if nargout==0, clear hl, end