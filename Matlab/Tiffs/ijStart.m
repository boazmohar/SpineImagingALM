function ijObj = ijStart
%ijStart: open ImageJ's main menu
%
% $Id: ijarray2plus.m 311 2008-07-31 20:25:12Z histed $

%global IMAGEJ_OBJ
%if isempty(IMAGEJ_OBJ)
%    ijObj = ij.ImageJ();
%    IMAGEJ_OBJ=ijObj;
%else
%    ijObj = IMAGEJ_OBJ;
%end
%ijObj = 
ijObj = ij.IJ.getInstance;
if isempty(ijObj)
    ijObj = ij.ImageJ();
end
ijObj.setTitle('IJ running in MATLAB JVM');
ijObj.show();
