function b=ijwritetiff(imp,fname)
%IJWRITETIFF Writes ImageJ ImagePlus object to disk
%   b=ijwritetiff(imp,fname)
%
if ~length(strfind(fname,'.tif'))
    fp=[fname,'.tif'];
else 
    fp = fname;
end

fs=ij.io.FileSaver(imp);
if imp.getImageStackSize==1
    b=fs.saveAsTiff(fp);
else    
    b=fs.saveAsTiffStack(fp);
end    

return;
