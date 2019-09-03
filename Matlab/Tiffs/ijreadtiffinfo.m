function [fi nImgInStack] = ijreadtiffinfo(pathname);
%IJREADTIFFINFO Reads TIFF file information and returns ImageJ TiffInfo object 
%   fi = ijreadtiffinfo(pathname);
%
% by Vincent Bonin based on original code by Aaron Kerlin
%
%   3/14/08 histed: support passing unix-style paths (forward slashes)
%
%$Id: ijreadtiffinfo.m 320 2008-08-22 13:58:05Z histed $

[directory,name,ext]=fileparts(pathname);

filename = [name ext];

if directory(end) ~= '\' || directory(end) ~= '/'
    directory(end+1) = filesep;
end
    
td=ij.io.TiffDecoder(directory,filename);
fiarray=td.getTiffInfo();


fi = fiarray(1);

% extract real number of images
nImgInStack = 0;
for iA=1:length(fiarray)
    nImgInStack = nImgInStack + fiarray(iA).nImages;
end
%  There seem to be two kinds of multipage tiffs.  One gives a list of
%  structures from getTiffInfo and one gives a single structure.  Here
%  nImgInStack handles both cases

return;
