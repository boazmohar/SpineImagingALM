function hl = fn_regression(x,y,varargin)
% function hl = fn_regression(x,y[,ha][,'square'])

% Thomas Deneux
% Copyright 2011-2012

% Input
ha = [];
dosquare = false;
for k=1:length(varargin)
    a = varargin{k};
    if ishandle(a)
        ha = a;
    elseif ischar(a)
        switch a
            case 'square'
                dosquare = true;
            otherwise
                error argument
        end
    else
        error argument
    end
end
if isempty(ha), ha = gca; end

% Fit
p = polyfit(x,y,1);
xx = [min(x) max(x)];
yfit = polyval(p,xx);

% Display
hl(1) = plot(x,y,'+','parent',ha);
hl(2) = line(xx,yfit,'color','k','linewidth',1.5,'parent',ha);
if nargout==0, clear hl, end

% Axis
if dosquare
    fn_axis(ha,'image',1.2)
else
    fn_axis(ha,'tight',1.2)
end
