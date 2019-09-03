function pos = fn_pixelpos(hobj)
% function pos = fn_pixelpos(hobj)
%---
% returns the position in pixels of any object without needing to
% change any units values
%
% See also fn_pixelsize 

% Thomas Deneux
% Copyright 2011-2012


switch get(hobj,'type')
    case 'figure'
        pos = get(hobj,'position');
    otherwise
        units = get(hobj,'units');
        switch units
            case 'pixels'
                pos = get(hobj,'position');
            case 'normalized'
                psiz = fn_pixelsize(get(hobj,'parent'));
                pos = get(hobj,'position').*psiz([1 2 1 2]);
            otherwise
                error('units ''%s'' not handled',units)
        end
end
 