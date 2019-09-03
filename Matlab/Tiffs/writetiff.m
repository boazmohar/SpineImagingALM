function b=writetiff(array,filename, typestr)
%WRITETIFF Writes array as TIFF file
%   b=writetiff(array,filename, typestr)
%
%   Uses ImageJ code to do write
%   typestr is a matlab type name, supplied to the cast function
%
%   array should be size [nRows,nCols,nFrames,nPlanes] where nPlanes is 1
%   (indexed color) or 3 (RGB color)
%
%   See also: readtiff
%
%$Id: writetiff.m 305 2008-07-18 20:55:10Z vincent $

tic;
if nargin < 3
    typestr = class(array);
end

if ~strcmp(class(array),typestr)
    array = cast(array,typestr);
end

[h,w,nframes,nPlanes] = size(array);

fprintf(1,'creating imageplus object\n');

imp = ijarray2plus(array,typestr);

fprintf(1,'writing file\n');

[pathstr]=fileparts(filename);

if length(pathstr) & exist(pathstr)~=7
    mkdir(pathstr);
end
    

b = ijwritetiff(imp,filename);

if ~b
    error('write fail: does directory exist?');
end

t= toc;
s=whos('array');
fprintf('wrote %i frames in %2.1f seconds (%2.0f fps or %5.0f MB/s )\n',nframes,t,nframes/t,s.bytes/1024^2/t);

clear imp; % probably not necessary


return;
