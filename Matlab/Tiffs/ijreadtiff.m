function imp = ijreadtiff(pathname);
%IJREADTIFF Reads TIFF file and returns ImageJ ImagePlus object
%   imp = ijreadtiff(pathname);
% 
%   Reads tif file and returns ImagePlus object
%
%   by Vincent Bonin based on original code by Aaron Kerlin
%
%   3/14/08 histed: support unix-style paths
%$Id: ijreadtiff.m 204 2008-04-30 19:07:26Z histed $ 

[directory,name,ext]=fileparts(pathname);

if directory(end) ~= '\' || directory(end) ~= '/'
    directory(end+1) = filesep;
end

td=ij.io.TiffDecoder(directory, [name ext]);
tfi=td.getTiffInfo();
op=ij.io.Opener();
imp=op.openTiffStack(tfi);

return;

