function ijObj = ijStop
%ijStop: stop all running instances of imageJ and clean up
%
%$Id: ijarray2plus.m 311 2008-07-31 20:25:12Z histed $

ijObj = ij.IJ.getInstance;
if isempty(ijObj)
    ijObj = ij.ImageJ();
else
    ijObj.quit()
end

