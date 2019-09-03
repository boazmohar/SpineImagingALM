function dendReadSWC(obj)
set(obj.handles.ustatH,'String','Loading swc')
[Name,Path] = uigetfile('*.SWC');
obj.cells.swcPath = [Path Name];
%%
[fid, ~] = fopen(obj.cells.swcPath);
if fid==-1
    throw(MException('DendriteSeg:swcCallback','Error opening file %s'))
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
Table = readtable(obj.cells.swcPath,'FileType','text','ReadVariableNames',0,'Delimiter',' ',...
    'HeaderLines',LineIgnore);
Table.Properties.VariableNames = VarNames;
%% up sample by 2
Table.x = round(Table.x);
Table.y = round(Table.y);
Table.z = round(Table.z+1); % z starts at 0 here
Table.r = Table.r;
obj.cells.dendriteTable = Table;