function imgOut = normImage(obj,imgIn)
MINMAX = obj.display.MINMAX;
brightset = obj.display.brightset;
imgIn(1,1) = MINMAX(1);
imgIn(1,2) = MINMAX(2);

imgIn(imgIn<MINMAX(1)) = MINMAX(1);
imgIn(imgIn>MINMAX(2)) = MINMAX(2);

% draw the new version
max_temp=max(max(imgIn));
imgOut=imgIn.*brightset;
imgOut(imgOut>max_temp) = max_temp;
