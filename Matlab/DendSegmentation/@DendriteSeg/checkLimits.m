function Error = checkLimits(obj,X,Y)
XLim = get(obj.handles.axH,'XLim');
YLim = get(obj.handles.axH,'YLim');
if X > XLim(1) && X < XLim(2) && Y > YLim(1) && Y < YLim(2)
    Error = false;
else
    Error = true;
end