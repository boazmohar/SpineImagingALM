%% convert Mat to Tiffs from Thunder
%% for means
Directory = '';
RunName = '';
cd(Directory)
List = dir([RunName '*.mat']);
List = {List.name};
for i = 1:length(List)
    Name = List{i};
    clear Mean;
    load(Name);
    outputFileName = [Name(1:end-3) 'tif'];
    if exist('Mean','var')
        writetiff(Mean,outputFileName);
    elseif exist('Target','var')        
        writetiff(Target,outputFileName);
    else
        disp('error')
    end
end
%%
clear all;
load('Run1expDict.mat')
% load('Run1ArrayOther.mat')
% CR = reshape(ArrayOther,expDict.cropDims(1), expDict.cropDims(2), expDict.cropDims(3)*size(ArrayOther,2));
% writetiff(CR,'Run1ArrayOther.tif')
load('Run1ArrayCorrectL.mat')
CR = reshape(ArrayCorrectL,expDict.cropDims(1), expDict.cropDims(2), expDict.cropDims(3)*size(ArrayCorrectL,2));
writetiff(CR,'Run1ArrayCorrectL.tif')
load('Run1ArrayCorrectR.mat')
CR = reshape(ArrayCorrectR,expDict.cropDims(1), expDict.cropDims(2), expDict.cropDims(3)*size(ArrayCorrectR,2));
writetiff(CR,'Run1ArrayCorrectR.tif')
cd Run1

load('Sp.mat')
ScanFields = length(Sp_s.All.x);
framerate = 1/(Sp_s.linesPerField/Sp_s.OptionsStruct.lineRate*ScanFields);
cd '..'
dlmwrite('Run1Dims.txt',double([expDict.cropDims framerate]), 'precision','%10.5f')