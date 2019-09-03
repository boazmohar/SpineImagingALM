function y = fn_timevector(x,t,outputtype)
% function times|count = fn_timevector(times|count,dt|tidx[,outputtype])
%---
% switch representation of a point process between the times of the events
% and the number of events per time bin
%
% Input:
% - count       vector of spike count at each time bin (or array or cell
%               array of vectors)
% - times       vector of spike times (or cell array thereof)
% - dt          scalar - size of time bin
% - tidx        vector of all time instants
% - outputtype  'times', 'rate' or 'count' [default: output type is the
%               opposite of input type]; use 'rateperperiod' or
%               'countperperiod' to count using time bins that are not
%               centered on the time points in tidx, but rather whose sides
%               are defined by the time points in tidx (in such case, the
%               number of bins will be one less than the number of time
%               points)
%
% Output:
% - count       column vector or array
% - times       row vector or cell array of row vectors

% Thomas Deneux
% Copyright 2010-2012

if nargin==0, help fn_timevector, return, end

% Multiple data?
xisarray = ~iscell(x) && all(size(x)>1);
if xisarray, x = squeeze(num2cell(x,1)); end % matrix -> cell array of column vectors
xismult = iscell(x);
if ~xismult, x = {x}; end
nx = numel(x);

% Input type
if xisarray
    inputtype = 'count';
else
    xtest = x{1};
    kfirstnz = find(xtest,1,'first');
    if isempty(kfirstnz) || mod(xtest(kfirstnz),1)
        inputtype = 'times';
    else
        inputtype = 'count';
    end
end

% Output type
if nargin<3
    outputtype = fn_switch(inputtype,'times','count','count','times');
elseif ~ismember(outputtype,{'times' 'count' 'rate' 'countperperiod' 'rateperperiod'})
    error('output type must be either ''times'', ''count'' or ''countperperiod''')
end
doperiod = any(strfind(outputtype,'perperiod'));
dorate   = any(strfind(outputtype,'rate'));
if ~strcmp(outputtype,'times'), outputtype = 'count'; end

% Conversion type
convtype = [inputtype(1) '2' outputtype(1)];
if doperiod && strcmp(inputtype,'count'), error('cannot count per period if input is already a count'), end

% Time specification
if isscalar(t)
    if doperiod
        error('when counting events inside defined periods, the time periods must be defined by a vector of time instants')
    end
    dt = t;
    t0 = 0;
    iseqspacing = true;
    if strcmp(inputtype,'times')
        nper = ceil((max([x{:}])-t0)/dt);
    else
        nper = length(x{1});
    end
    if strcmp(convtype,'c2t') && nper==0
        disp 'help me solve this case'
        keyboard
        error('times to count conversion: number of time instants unknown, please provide time information as the vector of time instants')
    end
else
    tidx = t;
    if doperiod
        % disconnected periods can be supplied as a (2 x nper) array of
        % time points
        isperiodsconnex = isvector(tidx);
        if ~isperiodsconnex
            if size(tidx,2)~=2, error('wrong format for periods'); end
            tidx = tidx'; tidx = tidx(:);
        end
    else
        if ~isvector(tidx), error('vector of time instants is not a vector!'), end
    end
    iseqspacing = all(abs(diff(tidx,2)<100*eps(class(tidx))));
    if ~doperiod && ~iseqspacing
        error('time points are not equally spaced')
    end
    nper = length(tidx) - doperiod;
    if iseqspacing
        dt = tidx(2)-tidx(1);
        if doperiod, t0 = tidx(1)+dt/2; else t0 = tidx(1); end
    end
end

% Prepare output
switch outputtype
    case 'time'
        y = cell(size(x));
    case 'count'
        y = zeros([nper size(x)]);
end

% Conversion
for k=1:nx
    xk = x{k};
    switch convtype
        case 't2t'
            y{k} = xk(:)';
        case 'c2c'
            y(:,k) = xk(:);
        case 't2c'
            times = xk;
            if iseqspacing
                % take advantage of the fact that time instants are equally spaced
                for i=1:length(times)
                    idx = 1+round((times(i)-t0)/dt);
                    if idx<=0 || idx>nper, continue, end
                    y(idx,k) = y(idx,k)+1;
                end
            else
                % find 'manually' to which period belongs each event
                for i=1:length(times)
                    idx = find(times(i)>=tidx,1,'last');
                    if isempty(idx) || idx>nper, continue, end
                    y(idx,k) = y(idx,k)+1;
                end                
            end
        case 'c2t'
            count = xk;
            times = zeros(1,sum(count));
            idx = 0;
            for i=1:length(count)
                ci = count(i);
                if ci
                    times(idx+(1:ci)) = t0+(i-1)*dt;
                    idx = idx+ci;
                end
            end
            y{k} = times;
    end
end

% Convert count to rate
if dorate
    if iseqspacing
        y = y / dt;
    else
        for k=1:nper
            y(k,:,:) = y(k,:,:) / (tidx(k+1)-tidx(k));
        end
    end
end

% Case of non-connex periods
if doperiod && ~isperiodsconnex, y = y(1:2:nper,:,:); end



