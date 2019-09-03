function [array,iframe] = readtiff(pathname,fileind,textInName,useLeicaNumbering,allowMissing,newType,stackLen);
%READTIFF Reads TIFF file or sequences of TIFF files into MATLAB array
%   ARRAY  = READTIFF(PATHNAME, FILEIND, TEXTINNAME, USELEICANUMBERING, ALLOWMISSING, newType);
%
%   assumes that files in alphabetical order and of same type
%
%   TEXTINNAME - if non-empty, read only the files that contain this
%       string, not case-sensitive.  
%       You can use additional regexp matching characters: help REGEXPI
%   USELEICANUMBERING - boolean, true if you are reading files output by
%       the Leica (odd numbering: Leica saves _t009_, then _t999_, 
%       then _t1000_ in file names)
%   ALLOWMISSING - boolean
%   NEWTYPE - string. cast stack to newType if not empty
%
%   [array,nframes] = readtiff(...) returns the number of frames loaded
%
% 08/01/20 can be called for stack and sequences
% 08/03/14 MH enable reading of Leica tiff seq files
% 08/10/30 VB returns number of frame loaded
% 25/03/10 VB can specify return type

%$Id: readtiff.m 637 2011-02-27 18:11:49Z vincent $

%% notes on leica:
% they use an odd naming convention: time series

%%% arg processing
if nargin < 2, fileind = []; end
if nargin < 3, textInName = ''; end
if nargin < 4, useLeicaNumbering = false; end
if nargin < 5, allowMissing = false; end
if nargin < 6, newType ='';end
if nargin < 7; stackLen = inf;end
    
%% get file list from path
if exist(pathname) == 2  % a single file 
    list = dir(pathname); % only used to get file sizes  
    selNs = 1;
    [pathstr,name,ext]=fileparts(pathname);
    filelist = {[name ext]};
    nFiles = 1;
    assert(nargin == 1, ...
           'First arg is a single file, pass only one parameter');
    sortNs=1;
    if isempty(pathstr)
        pathstr='.';
    end
    
elseif exist(pathname)  == 7 % directory
    %% list directory
    pathstr = pathname;
    list = dir(pathname); 
    list(1:2)=[];
    filelist = {list.name};
    if isempty(filelist)
        error('No files found at path %s', pathname);
    end
    nOrig = length(filelist);
    origNs = 1:nOrig;
    selNs = origNs;
    origList = filelist;

    %% first restrict by indices passed in
    if ~isempty(fileind)
        if fileind < 1, error('fileind must be >= 1 (1-origin numbering)'); end
        selNs = selNs(fileind);
        filelist = origList(selNs);
    end
    
    %% remove non-tiffs
    rMatchNs = regexpi(filelist, '\.(tif|tiff)$');
    selNs = selNs(~cellfun(@isempty, rMatchNs));
    filelist = origList(selNs);
    
    %% restrict list by text str if necessary
    if ~isempty(textInName)
        rMatchNs = regexpi(filelist, textInName);
        selNs = selNs(~cellfun(@isempty, rMatchNs));
        filelist = origList(selNs);
    end
    
    %--- have full list by the time we hit here

    % now sort
    nFiles = length(filelist);
    %% remap leica numbers if requested
    if useLeicaNumbering
        flRenum = regexprep(filelist, '^(.*_t)([0-9][0-9][0-9])(_ch.*)$', ...
                            '$10$2$3');
        [crap, sortNs] = sort(flRenum);
    else
        leicaFmtInd = regexp(filelist, 't[0-9]+_ch0[01]\.tif');
        isLeicaFmt = ~cellfun(@isempty,leicaFmtInd);
        if any(isLeicaFmt)
            error('Data looks like Leica files but useLeicaNs == false');
        end
        sortNs = 1:nFiles;
    end
    
    %% check number resulting
    if nFiles<1
        % stop here: either return an empty or raise an error
        if allowMissing
            array = []; 
            return
        else
            error('No files remaining after restrictions applied');
        end
    else
        fprintf(1,'Found %i files \n', nFiles);
    end

else  % could not find directory/file
    error('pathname does not exist: %s', pathname);
end

% make list of files with full path
for iF = 1:nFiles
    % full file list in order on disk
%     fullFileList{iF} = fullfileMH(pathstr, filelist{iF});  
    fullFileList{iF} = fullfile(pathstr, filelist{iF});    
end

%% Get number of total frames so we can pre-allocate memory
fprintf(1, 'Collecting frame information ');
% first get all sizes, to enable size optimization
allSizes = {list.bytes};
selSizes = cat(1,allSizes{selNs});
lastSize = NaN;
lastNImages = 0;

% read in order of sort above
try
    nframes = 0;
    for ifile = sortNs
        % check for same size; if so skip
        tSize = selSizes(ifile);
        if tSize == lastSize
            nImages = lastNImages;
        else
            % have to read it
            [fi nImages] = ijreadtiffinfo(fullFileList{ifile});   
            lastSize = tSize;
            lastNImages = nImages;
            fprintf('.');
        end
        nframes = nframes + nImages;

    end
catch
    le = lasterror;
    if strcmp(le.identifier, 'MATLAB:Java:GenericException')
        fprintf(1, ['\n\n*** %s: Java error reading file, case mismatch or ' ...
                    'OS-specific path used?\n\n'], ...
                mfilename, fullFileList{ifile});
        ds = dir(fullFileList{ifile}) % disp directory entry
        if ds.bytes == 0
            error('TIFF file is zero bytes\n');
        end
    end
    rethrow(lasterror);
end


    

fprintf('\nBit depth  = %d. Number of frames = %d\n',fi.getBytesPerPixel*8, nframes);

% not sure what this is supposed to achieve 3/14/08 - MH
%  (note that fnl always has length 1 given the loop above)
% $$$ if length(unique(fnl)) > 1
% $$$     disp(fn);
% $$$     error('files not in alphabetical order');
% $$$ end

switch fi.fileType
    case fi.GRAY8
        dataType = 'uint8';
    case fi.GRAY16_SIGNED
        dataType = 'uint16';
    case fi.GRAY16_UNSIGNED
        dataType = 'uint16';
    case fi.GRAY32_FLOAT
        dataType = 'single';
    case fi.GRAY32_UNSIGNED
        dataType = 'uint32';
    otherwise
        error('unknown tiff file type');
end

if isempty(newType);    newType = dataType;end

if ~isfinite(stackLen); stackLen = nframes;end;

% pre-allocate stack of size stackLen
array = zeros(fi.height,fi.width,stackLen,newType);

tic;
totalT = toc;
elT = toc;
iframe = 0;
fprintf(1, 'Reading %d tiff files from %s\nframe ', nFiles, pathstr);
for ifile  = 1:nFiles
    % read in sorted order
    tFrN = sortNs(ifile);
    %disp(filelist{tFrN});
    imp = ijreadtiff(fullFileList{tFrN});        
    r1 = ijplus2array(imp,newType);
    
    if stackLen >= iframe+imp.getImageStackSize
        array(:,:,iframe+1:iframe+imp.getImageStackSize) = r1;
    else
        excess = (iframe + imp.getImageStackSize) - stackLen;
        array(:,:,iframe+1:stackLen) = r1(:,:,1:imp.getImageStackSize-excess);
        iframe = iframe + imp.getImageStackSize-excess;
        break;
    end
    
    iframe = iframe +  imp.getImageStackSize;        
    clear imp;    

    tElapsedS = toc - elT;
    if tElapsedS > 2,
        fprintf(1, '%d ', iframe); % t
        elT = toc;  
    end
end


s=whos('array');

totS = toc-totalT;
fprintf('\nloaded %i frames in %2.1f seconds (%2.0f fps or %5.0f MB/s )\n', ...
        iframe, totS, iframe/totS,s.bytes/1024^2/totS);
        
return;


% $$$ %% test readtiff
% $$$ pathname = 'I:\users\vincent\images\cat080116\cat080116A1_\cat080116A1__green';
% $$$ array = readtiff(pathname);
% $$$ pathname = 'I:\users\vincent\analysis\cat080116B1_\Registered.tif';
% $$$ array = readtiff(pathname);
