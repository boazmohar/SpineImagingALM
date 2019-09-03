function h = fn_hash(inp,varargin)
% function h = fn_hash(inp[,meth][,n])
%---
% This function relies on function hash.m downloaded on Matlab File
% Exchange.
% It extends hash functionalities to structures and cell arrays
% default method 'MD2' is used, input 'n' outputs n letters
% 
% structure with fields in different order give the same hash key

% Michael Kleder (function hash.m)
% Copyright 2005-2012
% Thomas Deneux
% Copyright 2007-2012

meth = 'MD2'; n = 0;
for i=1:length(varargin)
    a = varargin{i};
    if ischar(a)
        meth = a;
    else
        n = a;
    end
end

if isempty(inp), inp = class(inp); end
if isnumeric(inp) || islogical(inp) || ischar(inp)
    h = hash(inp,meth);
else
    % transform object into cell array with numeric / character type elements
    if isstruct(inp) || isobject(inp)
        inp = orderfields(struct(inp));
        F = fieldnames(inp);
        C = struct2cell(inp);
        C = [F(:); C(:)];
    elseif iscell(inp)
        C = inp(:);
    else
        error('cannot hash object of class ''%s''',class(inp))
    end
    C = cat(1,C,{class(inp); size(inp)});
    
    % hash each element of the cell array
    h = cell(1,numel(C));
    for i=1:numel(C)
        h{i} = fn_hash(C{i},meth);
    end
    
    % hash the concatenation of all hash results
    h = hash([h{:}],meth);
end

% make a n-letters word out of it
if n
    h = h(1:n);
    f = (h>='0' & h<='9');
    h(f) = h(f)-'0'+'A';
    h(~f) = h(~f)-'a'+'K';
end

%---
function h = hash(inp,meth)
% HASH - Convert an input variable into a message digest using any of
%        several common hash algorithms
%
% USAGE: h = hash(inp,'meth')
%
% inp  = input variable, of any of the following classes:
%        char, uint8, logical, double, single, int8, uint8,
%        int16, uint16, int32, uint32, int64, uint64
% h    = hash digest output, in hexadecimal notation
% meth = hash algorithm, which is one of the following:
%        MD2, MD5, SHA-1, SHA-256, SHA-384, or SHA-512 
%
% NOTES: (1) If the input is a string or uint8 variable, it is hashed
%            as usual for a byte stream. Other classes are converted into
%            their byte-stream values. In other words, the hash of the
%            following will be identical:
%                     'abc'
%                     uint8('abc')
%                     char([97 98 99])
%            The hash of the follwing will be different from the above,
%            because class "double" uses eight byte elements:
%                     double('abc')
%                     [97 98 99]
%            You can avoid this issue by making sure that your inputs
%            are strings or uint8 arrays.
%        (2) The name of the hash algorithm may be specified in lowercase
%            and/or without the hyphen, if desired. For example,
%            h=hash('my text to hash','sha256');
%        (3) Carefully tested, but no warranty. Use at your own risk.
%        (4) Michael Kleder, Nov 2005
%
% EXAMPLE:
%
% algs={'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
% for n=1:6
%     h=hash('my sample text',algs{n});
%     disp([algs{n} ' (' num2str(length(h)*4) ' bits):'])
%     disp(h)
% end

inp=inp(:);
% convert strings and logicals into uint8 format
if ischar(inp) || islogical(inp)
    inp=uint8(inp);
else % convert everything else into uint8 format without loss of data
    inp=typecast(inp,'uint8');
end

% verify hash method, with some syntactical forgiveness:
meth=upper(meth);
switch meth
    case 'SHA1'
        meth='SHA-1';
    case 'SHA256'
        meth='SHA-256';
    case 'SHA384'
        meth='SHA-384';
    case 'SHA512'
        meth='SHA-512';
    otherwise
end
algs={'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
if isempty(strmatch(meth,algs,'exact'))
    error(['Hash algorithm must be ' ...
        'MD2, MD5, SHA-1, SHA-256, SHA-384, or SHA-512']);
end

% create hash
x=java.security.MessageDigest.getInstance(meth);
x.update(inp);
h=typecast(x.digest,'uint8');
h=dec2hex(h)';
if(size(h,1))==1 % remote possibility: all hash bytes < 128, so pad:
    h=[repmat('0',[1 size(h,2)]);h];
end
h=lower(h(:)');
clear x
return
