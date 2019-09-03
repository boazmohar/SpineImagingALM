function [bwNew, bwCell, cellMean] = subAddCell(obj, bwCurr, fovImg, cellMean)
%% add pixels in the current Z plane
X = obj.current.X;
Y = obj.current.Y;
cellRadius = 1;
diskR= obj.current.diskR;
cThresh = obj.current.cThresh;
[nRows, nCols] = size(bwCurr);
%% set up dilation disks:  use 4-connected regions for all
se = obj.display.se;  % for cells-to-avoid: all 9 pix in a square;
                       % avoid diagonally-connected cells.
se2 = obj.display.se2;   % for region to find mean over
seJunk = strel('disk', max(round(cellRadius/4), 1), 4);  % remove thin junk
seExpand = strel('disk', diskR, 4);  % expand thresholded region
%% add a disk around each point, non-overlapping with adj cells
tempmask = false(nRows, nCols);
dilateorg = imdilate(bwCurr,se);
tempmask(Y, X) = 1;
tempmask = imdilate(tempmask,se2);
tempmask = tempmask & ~dilateorg;
%% fill region around disk of similar intensity, combine with disk
if nargin < 4
    cellMean = mean(fovImg(tempmask),1);
end
%%
if nargout > 1
    if fovImg(Y, X) > cellMean.*cThresh
        allMeanBw = fovImg >= cellMean.*cThresh;  % threshold by intensity
        %Npts = sum(sum(allMeanBw))
        connMeanBw = bwselect(allMeanBw &~dilateorg, X, Y, 4);
        connMeanBw = connMeanBw |tempmask & ~dilateorg;

        % erode then dilate filled to remove sharp things
        erMean = imerode(connMeanBw, seJunk);
        dilateMean = imdilate(erMean, seJunk);
        dilateMean = imdilate(dilateMean, seExpand); % (thresh is conservative)
        bwCell = dilateMean & ~dilateorg;
        bwNew = bwCurr | bwCell;
    else
        bwCell = false(size(bwCurr));
        bwNew = bwCurr;
    end
end
