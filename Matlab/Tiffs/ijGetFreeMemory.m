function ijGetFreeMemory(object, event)
%IJGETFREEMEMORY prints ImageJ Java Heap usage to standard input
%   ijGetFreeMemory(object, event)

fprintf(1,'Java Heap Usage %s\n',char(ij.IJ.freeMemory));

return
