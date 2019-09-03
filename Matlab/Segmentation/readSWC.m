function Table = readSWC(Options)
% Table = readSWT(Options) retunes a table from a SWC file
% Options is ia strucy contanint the filename to read and 
% typeNumberss which is an array of types that will be returned from the file
% if typeNumberss = -1 all types will be returned
% if typeNumberss is ommited defulat is 1
if nargin < 1
    throw(MException('MATLAB:ambiguousSyntax','Specify Options'))
end
if isfield(Options,'filename')
    filename = Options.filename;
else
    throw(MException('MATLAB:ambiguousSyntax','Specify filename'))
end
if isfield(Options,'typeNumbers')
    typeNumbers = Options.typeNumbers;
else
    typeNumbers = 1;
end
[fid, ~] = fopen(filename);
if fid==-1
    throw(MException('MATLAB:ambiguousSyntax','Error opening file'))
end
tline = fgets(fid);
LineIgnore =0;
while ischar(tline) && tline(1) == '#'
    LineIgnore = LineIgnore+1;
    tline = fgets(fid);
end
fclose(fid);
VarNames = {'Num','Type','x','y','z','r','pNum'};
Table = readtable(filename,'FileType','text','ReadVariableNames',0,'Delimiter',' ',...
    'HeaderLines',LineIgnore);
Table.Properties.VariableNames = VarNames;
if typeNumbers ~= -1
    for i = 1:length(typeNumbers)
        if i == 1
            Rows = Table.Type == typeNumbers(i);
        else
            Rows = Rows | Table.Type == typeNumbers(i);
        end
    end
    if sum(Rows) > 0
        Table = Table(Rows,:); 
    else
        disp(['No type: ' num2str(typeNumbers) ' in the file: ' filename]);
    end
end