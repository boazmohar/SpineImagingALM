function [opts,opts2]=optndfts(args,varargin)
% optndfts  Create option struct with default values for those not given
%  2015-09-07  Matlab  Copyright (c) 2015, W J Whiten  BSD License
%
% opts=optndfts(args,dflts)
% opts=optndfts(args,initial,dflts)
%  args     Typically varargin from calling function
%            A cell array containing: initial arguments if present  
%            then name value pairs, or containing a struct,
%            or containing a cell array of name value pairs.
%            Alternatively a struct. Note initial arguments can not be of
%            type char, and a single struct or cell, in a cell array 
%            will be interpreted as giving option names and values.
%  initial  A cell array of the names of optional initial arguments given
%            without names in varargin of the calling function (ie args).
%            These names are the given order for the optional initial 
%            arguments and zero or more of the initial arguments can
%            can be given in args, (the optional initial arguments 
%            in args can not be of type char (unless preceeded by -), or
%            a be single struct, single cell array or single empty array,
%            these cases must be given as a name value pair).
%            The names given here must have default values in dflts.
%            This argument is optional.
%  dflts    Argument list giving default name value pairs, or a struct, 
%            or a cell array of name value pairs.
%
%  opts     Struct giving names and values for all variables in dflts
%  opts2    Struct giving names & values of fields not in opts (optional)
%            If given fields not in dflts do not give an error
%            Can be used to pass options to nested functions (note
%            alternative is to put in a field containing a cell or struct 
%            with the options for nested functions)
%
% Use to process optional arguments given as name value pairs, a cell 
%  array of name value pairs, or as a struct. Default values from 
%  optndfts call are used for values not given in args.
% Note avoid struct arrays eg struct('a',{1 2}); use struct('a',{{1 2}})
% 
% Example:
% Using optndft within a function as follows: 
%    function strt=test1(a,varargin) % varargin values replace defaults
%    strt=optndfts(varargin,'aa',1,'bb',2,'cc',3); % give default values
%    end
% Then demo could be called as:
%    test1(123,'bb',22)
% And optndfts fills in default values for those not given::
%    ans = 
%        aa: 1
%        bb: 22
%        cc: 3
% 
% Also can allow for optional arguments that are not named:
%    function strt=test2(arg1,varargin)  % varargin replaces defaults
%    strt=optndfts(varargin,{'aa','bb'},'aa',1,'bb',2,'cc',3);
%    return  	% optndfts returns a struct
%    end
% 
%    test2(123,11,'cc',33)
%    ans = 
%        aa: 11
%        bb: 2
%        cc: 33
% 
% See also:  inputParser

% check for two or three arguments
if(length(varargin)>1 && iscellstr(varargin{1}))
    initial=varargin{1};
    varargin=varargin(2:end);
else
    initial={};
end

% convert varargin inputs to a struct
if(length(varargin)==1 && iscell(varargin))
    varargin=varargin{1};
end
if(iscell(varargin))
    opts=struct();
    for i=1:2:length(varargin)
        opts.(varargin{i})=varargin{i+1};
    end
else
    opts=varargin;
end

% initial value for second output
opts2=struct();

% if input args not a struct convert to struct
firstargs={};
if(iscell(args))
    
    % first check for initial arguments without name
    i=0;
    if(~(length(args)==1 &&  ...
            (isstruct(args{1}) || iscell(args{1}) || isempty(args{1}))))
        while(i<length(args))
            if(ischar(args{i+1}))
                if(args{i+1}(1)=='-')
                    args{i+1}=args{i+1}(2:end);
                else
                    break
                end
            end
            i=i+1;
        end
    end
    if(i>length(initial))
        error('optndfts: Too many unnamed input arguments')
    end
    firstargs=args(1:i);
    args=args(i+1:end);
    
    if(length(args)==1 && iscell(args))
        args=args{:};
    end
    if(iscell(args))    % for cell array of name value pairs
        args1=struct();
        for i=1:2:length(args)
            args1.(args{i})=args{i+1};
        end
        args=args1;
    end
end

% initial arguments without names in args, copy over defaults
for i=1:length(firstargs)
    if(~isfield(opts,initial{i}))
        error(['optndfts: Initial name  ',initial{i},  ...
            '  does not have default value'])
    else
        opts.(initial{i})=firstargs{i};
    end
end

% copy new values over default values
if(isempty(args))
    args=struct();
end
names=fieldnames(args);
for i=1:length(names)
    if(~isfield(opts,names{i}))
        if(nargout>1)
            opts2.(names{i})=args.(names{i});
        else
            error(['optndfts: Invalid field name  ',names{i}])
        end
    end
    opts.(names{i})=args.(names{i});
end

return
end
