function s = fn_readxml(fname)
% function s = fn_readxml(fname)
%---
% This function is a wrapper for function xml_read.m downloaded on Matlab
% File Exchange.
% Additionally, the structure is slightly simplified.

% Jarek Tuszynski (xml_read.m function)
% Copyright 2007-2012
% Thomas Deneux
% Copyright 2007-2012

if nargin==0, help fn_readxml, return, end

s = xml_read(fname);
s = simplify(s);

%--- recursive sub-function ---%
function s = simplify(s)

F = fieldnames(s);
for k=1:length(F)
    f = F{k};
    if strcmp(f,'ATTRIBUTE')
        error 'please check if it makes more sense to alter the structure inside fn_readxml, or inside the calling function'
        F2 = fieldnames(s(1).ATTRIBUTE);
        for k2=1:length(F2)
            f2 = F2{k2};
            for i=1:numel(s)
                s(i).(f2) = s(i).ATTRIBUTE.(f2);
            end
        end
        s = rmfield(s,'ATTRIBUTE');
    elseif isstruct(s(1).(f))
        for i=1:numel(s)
            if ~isempty(s(i).(f))
                s(i).(f) = simplify(s(i).(f));
            end
        end
    end
end


