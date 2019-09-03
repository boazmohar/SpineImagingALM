function x = fn_minmax(varargin)
% function fn_minmax(action,x1,x2,..)
%
% utility to easily compute best axes positions...
% 'action' determines what to do :
% - 'minmax'    -> [min(xi(1)) max(xi(2))]
% - 'axis'      -> [min max min max]
% - logical vector -> ok

% Thomas Deneux
% Copyright 2006-2012

if nargin==0, help fn_minmax, return, end

action = varargin{1};
if ischar(action)
    switch action
    case 'minmax'
        action = [0 1];
    case 'axis'
        action = [0 1 0 1];
    end
end
if ~all(action==0 | action==1), error('bad action argument'); end

if any(size(action)==1)
    x = zeros(length(action),1)*nan;
    imax = find(action);
    imin = find(~action);
    if any(size(varargin{2})==1)
        for i=1:nargin-1 % vectors
            xi = varargin{i+1}(:);
            if ~any(size(xi)==1) | length(xi)~=length(action)
                error('arguments are not the same length')
            end
            x(imin) = min(x(imin),xi(imin));
            x(imax) = max(x(imax),xi(imax));
        end
    else % one matrix
        x1 = varargin{2};
        if nargin>2, error('too many arguments');
            x(imin) = min(x1(:,imin));
            x(imax) = max(x1(:,imax));
        end
    end
else % matrices
    x = zeros(size(action))*nan;
    imax = find(action);
    imin = find(~action);
    for i=1:nargin-1
        xi = varargin{i+1};
        if size(xi)~=size(action)
            error('arguments are not the size')
        end
        x(imin) = min(x(imin),xi(imin));
        x(imax) = max(x(imax),xi(imax));
    end
end
