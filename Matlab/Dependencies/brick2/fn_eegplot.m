function hl = fn_eegplot(varargin)
% function hl = fn_eegplot(usual plot arguments[,stepflag]['num'])
%---
% Like Matlab function 'plot', but separates line by a distance 'ystep'
% the way 'ystep' is calculated is determinated by 'stepflag':
%   - x         use ystep = x
%   - 'STD'     use ystep = mean(std(data))
%   - 'xSTD'    use ystep = x*mean(std(data)) [default is STD]
%   - 'fit' 	use ystep = max(max(data)-min(data))
%   - 'xfit'    use ystep = x * max(max(data)-min(data))
%
% additional flag 'num' indicates to rescale the data so that numbers on
% the ordinates axis correspond to data number

% Thomas Deneux
% Copyright 2005-2012

% Input
% (data at position 1 or 2)
if nargin>1 && isnumeric(varargin{2})
    idata = 2;
else
    idata = 1;
end
data = varargin{idata};

% (flags at the end)
donum = false;
ystep = [];
% ('num' flag?)
a = varargin{end};
if ischar(a) && strcmp(a,'num')
    donum = true; 
    varargin(end) = [];
end
% (step specification)
a = varargin{end};
if isnumeric(a) && length(a)==1 && ~(ishandle(a) && strcmp(get(a,'type'),'axes')) % TODO: not enough!!
    ystep = a;
    varargin(end) = [];
elseif ischar(a)
    x = regexpi(a,'^([0-9\.]*)STD$','tokens');
    if ~isempty(x)
        if isempty(x{1}{1}), fact=1; else fact=str2double(x{1}); end
        ystep = fact*mean(std(data));
        varargin(end) = [];
    else
        x = regexp(a,'^([0-9\.]*)fit$','tokens');
        if ~isempty(x)
            if isempty(x{1}), fact=1; else fact=str2double(x{1}); end
            ystep = fact*max(max(data(:))-min(data(:)));
            varargin(end) = [];
        end
    end
end

% more computation
if isempty(ystep)
    ystep = mean(std(data));
end
if donum
    data = data/ystep;
    varargin{idata} = data;
    ystep = 1;
end

% display
hl = plot(varargin{:});
for k=1:length(hl)
    set(hl(k),'ydata',get(hl(k),'ydata')+(k-1)*ystep)
end
axis tight

if nargout==0, clear hl, end

    
    
    
    