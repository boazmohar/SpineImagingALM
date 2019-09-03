function s = fn_structmerge(s,s1,varargin)
% function s = fn_structmerge(s,s1[,'skip|strict'][,'recursive'][,'type'])
%---
% set or replace values in s from those in s1, where s and s1 are
% structures of the same size
% - if 'skip', or 'strict' flag is specified: does not add new field in
%   structure s (generates error if 'strict' flag and s1 has additional
%   fields)
% - if 'recursive' flag, field values which are themselves structures are
%   not merely replaced, but are also merged using fn_structmerge 
% - if 'type' flag is specified: also requires field values to be the same
%   class in s and s1 when the field already exists in s (generates error
%   if it is not the case, except that it performs the conversions
%   0/1->false/true  and char array->cell array of strings)

% Thomas Deneux
% Copyright 2007-2012

if isempty(s1), return, end
if ~isscalar(s1) && any(size(s1)~=size(s))
    error('size mismatch')
end
[skip strict recursive type] = fn_flags('skip','strict','recursive','type',varargin);
skip = skip | strict;

F = fieldnames(s1);
for k=1:length(F)
    f = F{k};
    if skip && ~isfield(s,f)
        if strict
            error('field ''%s'' not present in original structure',f)
        end
    elseif recursive && isfield(s,f) && isstruct(s(1).(f))
        for i=1:numel(s)
            if isscalar(s1), j=1; else j=i; end
            s(i).(f) = fn_structmerge(s(i).(f),s1(j).(f),varargin{:});
        end
    else
        for i=1:numel(s)
            if isscalar(s1), j=1; else j=i; end
            if type && isfield(s,f) && ~strcmp(class(s(i).(f)),class(s1(j).(f)))
                if islogical(s(i).(f)) && isnumeric(s1(j).(f)) ...
                        && isscalar(s1(j).(f)) && ismember(s1(j).(f),[0 1])
                    s(i).(f) = logical(s1(j).(f));
                elseif iscell(s(i).(f)) && (isempty(s(i).(f)) || ischar(s(i).(f){1})) ...
                        && ischar(s1(i).(f))
                    s(i).(f) = cellstr(s1(i).(f));
                else
                    error('class mismatch')
                end
            else
                s(i).(f) = s1(j).(f);
            end
        end
    end
end