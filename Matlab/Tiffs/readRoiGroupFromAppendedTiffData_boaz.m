function [hMroiRoiGroup hStimRoiGroups] = readRoiGroupFromAppendedTiffData_boaz(filename)
    hMroiRoiGroup = [];

    if nargin ~= 1 
        disp('Expected filename as argument');
        return;
    end

    [roiStr numBytes] = readAppendedString(filename);
	appendedData = scanimage.mroi.util.deserialize(roiStr, 1);

	%MROI group deserialization
	if isfield(appendedData,'mroiRoiGroup') 
		hMroiRoiGroup = appendedData.mroiRoiGroup;
	end

	%Photostim groups deserialization
	if isfield(appendedData,'stimRoiGroups') 
		if numel(appendedData.stimRoiGroups) == 0
			hStimRoiGroups = [];
		else
			hStimRoiGroups = appendedData.stimRoiGroups;
		end
	else
		hStimRoiGroups = [];
	end
end

function [str, numBytes] = readAppendedString(filename)
    %readAppendedString - reads non-TIFF data appended string from a file
    % Usage :	writeraw(G, filename)
    % str:		string to append to the image
    % filename: file name of the file to append to 
    % count:	return value, the elements written to file

    str = [];
    %disp(['Reading '  filename ' ...']);

    % Get file ID
    fid = fopen(filename,'r');

    % Check if file exists
    if (fid == -1)
        error('Cannot open file\n');
        pause
    end

    %Move to the end of file
    numBytesSize = 4;
    ret = fseek(fid, -numBytesSize, 1);
    if ret ~= 0 
        disp('Error: fseek failed');
        fclose(fid);
        return;
    end
    numBytes = fread(fid,numBytesSize,'uint32');

    %Move to the beginning of the desired chunk
    ret = fseek(fid, -numBytesSize-numBytes, 1);
    if ret ~= 0 
        disp('Error: fseek failed');
        fclose(fid);
        return;
    end
    %str = fread(fid,numBytes,'uint8=>char');
    str = fread(fid,numBytes,'*char')';


    fclose(fid);
    %disp(['Closing '  filename ' ...']);

end %function

