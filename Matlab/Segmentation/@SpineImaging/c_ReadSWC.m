function Table = c_ReadSWC(obj)
% Table = readSWT(obj) retunes a table from a SWC file
%% open file to see how many lines to skip
filename = obj.OptionsStruct.filename;
[fid, ~] = fopen(filename);
if fid==-1
    throw(MException('SpineImaging:readSWC','Error opening file %s'))
end
tline = fgets(fid);
LineIgnore =0;
while ischar(tline) && tline(1) == '#'
    LineIgnore = LineIgnore+1;
    tline = fgets(fid);
end
fclose(fid);
%% read the table
VarNames = {'Num','Type','x','y','z','r','pNum'};
Table = readtable(filename,'FileType','text','ReadVariableNames',0,'Delimiter',' ',...
    'HeaderLines',LineIgnore);
Table.Properties.VariableNames = VarNames;
%% Converting pixels to um
UM_Zoom = obj.OptionsStruct.UM_1X./obj.OptionsStruct.Zoom;
Pixel_Size = UM_Zoom./obj.OptionsStruct.RefPixels;
Table.x = Table.x.*Pixel_Size;
Table.y = Table.y.*Pixel_Size;
Table.z = Table.z.*obj.OptionsStruct.Zstep;
%% setting the object's table
obj.Table = Table;