% [cellMask cellBoundary] = getCellMaskMSS(img, varargin) estimates the cell mask/outline using multi-scale steerable filters
%
% Inputs:
%             img : input image
% 
% Options:
%        'Scales' : vector of scales (sigma) used by the filter. Default: [1 2 4].
%   'FilterOrder' : order of the filters. Default: 3.
%  'RemoveRadius' : radius of the final erosion/refinement step
%
% Outputs:
%        cellMask : binary mask of the cell 
%    cellBoundary : binary mask of the cell outline
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
%
% This file is part of WindowingPackage.
% 
% WindowingPackage is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% WindowingPackage is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with WindowingPackage.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Francois Aguet, September 2011 (last modified: 10/23/2011)

function [cellMask, cellBoundary] = getCellMaskMSS(img, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addParamValue('Scales', [1 2 4], @isvector);
ip.addParamValue('FilterOrder', 3, @(x) ismember(x, [1 3 5]));
ip.addParamValue('SearchRadius', 6, @isscalar);
ip.addParamValue('NormalizeResponse', false, @islogical);
ip.addParamValue('Mask', []);
ip.parse(img, varargin{:});
scales = ip.Results.Scales;

[ny,nx] = size(img);
% ordered index, column-order CCW
borderIdx = [1:ny 2*ny:ny:(nx-1)*ny nx*ny:-1:(nx-1)*ny+1 (nx-2)*ny+1:-ny:ny+1];
borderMask = zeros(ny,nx);
borderMask(borderIdx) = 1;

%------------------------------------------------------------------------------
% I. Multi-scale steerable filter
%------------------------------------------------------------------------------
[res, theta, nms, scaleMap] = multiscaleSteerableDetector(img, ip.Results.FilterOrder, scales);

if ip.Results.NormalizeResponse
    res = res ./ filterGauss2D(res, 5);
    %nms = nonMaximumSuppression(res, theta);
    % -or-
    nms = (nms~=0).*res; % better
end

% Mask of candidate edges
edgeMask = double(bwmorph(nms~=0, 'thin'));

% Break any Y or higher order junctions
nn = (imfilter(edgeMask, ones(3), 'same')-1) .* edgeMask;
junctionMask = nn>2;
segmentMask = edgeMask .* ~junctionMask;

% endpoints of all segments
% endpointMatrix = (imfilter(segmentMatrix, ones(3), 'same')-1) .* segmentMatrix;
% endpointMatrix = endpointMatrix==1;

% generate list of segments and add associated properties
CC = bwconncomp(segmentMask, 8);

% identify and remove single pixels, update segment matrix
csize = cellfun(@numel, CC.PixelIdxList);
singletonIdx = CC.PixelIdxList(csize==1);
segmentMask([singletonIdx{:}]) = 0;
CC.PixelIdxList(csize==1) = [];
CC.NumObjects = numel(CC.PixelIdxList);
csize(csize==1) = [];
CC.NumPixels = csize;

% labels of connected components
labels = double(labelmatrix(CC));

% compute pixel order and intensity on each side for all detected edge segments
CC = computeEdgeSegmentProperties(CC, img, theta);

%hval = zeros(1,CC.NumObjects);
%for k = 1:CC.NumObjects
%    hval(k) = kstest2(CC.rval{k}(:), CC.lval{k}(:));
%end

%------------------------------------------------------------------------------
% II. Rough estimate of the cell outline based on threshold: coarseMask
%------------------------------------------------------------------------------
coarseMask = ip.Results.Mask;
if isempty(coarseMask)
    % threshold 1st mode (background) of histogram
    img_smooth = filterGauss2D(img, 1);
    T = thresholdFluorescenceImage(img_smooth);
    coarseMask = double(img_smooth>T);
    coarseMask = bwmorph(coarseMask, 'fill'); % clean up isolated negative pixels
end
% get boundary from this mask
bdrVect = bwboundaries(coarseMask);
bdrVect = vertcat(bdrVect{:});
coarseBdr = zeros(ny,nx);
coarseBdr(sub2ind([ny nx], bdrVect(:,1), bdrVect(:,2))) = 1;

% endpoints/intersection of boundary w/ border
borderIS = coarseBdr & borderMask;
borderIS = double(borderIS(borderIdx));
borderIS = borderIdx((conv([borderIS(end) borderIS borderIS(1)], [1 1 1], 'valid')-1)==1);

% clean up, remove image border, add back border intersections
coarseBdr = bwmorph(coarseBdr, 'thin');
coarseBdr(borderIdx) = 0;
coarseBdr(borderIS) = 1;

edgeSearchMask = imdilate(coarseBdr, strel('disk', 20)); % ! dilation is arbitrary !

% labels within search area
idx = unique(labels.*edgeSearchMask);
idx(idx==0) = []; % remove background label

% update connected components list
CC.NumObjects = numel(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
CC.NumPixels = CC.NumPixels(idx);
CC.PixelOrder = CC.PixelOrder(idx);
CC.EndpointIdx = CC.EndpointIdx(idx);
CC.rval = CC.rval(idx);
CC.lval = CC.lval(idx);
% labels = double(labelmatrix(CC));

% order interpolated intensities such that the lower intensities are always in 'lval'
for k = 1:CC.NumObjects
    rmean = nanmean(CC.rval{k}(:));
    lmean = nanmean(CC.lval{k}(:));
    if rmean<lmean
        tmp = CC.rval{k};
        CC.rval{k} = CC.lval{k};
        CC.lval{k} = tmp;
    end
end

% mask with average intensity of each segment
avgInt = cellfun(@(px) sum(nms(px)), CC.PixelIdxList) ./ CC.NumPixels;
edgeMask = zeros(ny,nx);
for k = 1:CC.NumObjects
    edgeMask(CC.PixelIdxList{k}) = avgInt(k);
end

% Some of the edges are background noise -> bimodal distribution of edge intensities
val = nms(edgeMask~=0); % intensities of edges
minv = min(val);
maxv = max(val);
T = graythresh(scaleContrast(val, [], [0 1]));
T = T*(maxv-minv)+minv;

% initial estimate of cell contour
cellBoundary = edgeMask > T;
%cellBoundary = edgeMask > min(T, thresholdRosin(val));
%cellBoundary = edgeMask;

% 1st graph matching based on orientation at endpoints, with small search radius
[matchedMask] = matchSegmentEndPoints(cellBoundary, theta, 'SearchRadius', ip.Results.SearchRadius, 'Display', false);
matchedMask = double(matchedMask);

% merge edge information: intensities on both sides, endpoints
CC = updateEdgeInfo(matchedMask, CC);

% The connected components in this mask are no longer simple segments. For the next matching 
% steps, the two outermost endpoints are needed.
CC.EndpointIdx = getEndpointsGeodesic(matchedMask, horzcat(CC.EndpointIdx{:})');

%avgInt = cellfun(@(px) sum(res(px)), CC.PixelIdxList) ./ CC.NumPixels;
% for k = 1:CC.NumObjects
%     %matchedMask(CC.PixelIdxList{k}) = avgInt(k);
% end

%----------------------------------------
% Matching based on edge similarity
%----------------------------------------
matchesFound = true;
iter = 0;
while matchesFound
    [matchLabel, matchIndex] = matchEdgeSegments(CC, 'SearchRadius', 10);
    if isempty(matchIndex)
        matchesFound = false;
    end
    matchedMask = double(linkEdgeSegments(matchedMask, matchIndex));
    CC = updateEdgeInfo(matchedMask, CC);
    CC.EndpointIdx = getEndpointsGeodesic(matchedMask, horzcat(CC.EndpointIdx{:})');
    %matchLabel
 
    iter = iter + 1;
end


%----------------------------------------
% Matching based on edge strength
%----------------------------------------
% [backMask,backMu,backSig] = estimateBackgroundArea(img);
% bgDist = img(backMask==1);

% check for connected components, whether further matching is required
% tmp = labels==0;
% CCX = bwconncomp(tmp,4)

matchesFound = true;
iter = 0;
%while matchesFound
    [matchLabel, matchIndex] = matchEdgeSegments(CC, 'SearchRadius', 10, 'Mode', 'Edge');
    if isempty(matchIndex)
        matchesFound = false;
    end
    matchedMask = double(linkEdgeSegments(matchedMask, matchIndex));
    CC = updateEdgeInfo(matchedMask, CC);
    CC.EndpointIdx = getEndpointsGeodesic(matchedMask, horzcat(CC.EndpointIdx{:})');
    %matchLabel

    iter = iter + 1;
%end



img0 = scaleContrast(img);
img1 = img0;
img0(matchedMask~=0) = 0;
img1(matchedMask~=0) = 255;
rgb = uint8(cat(3, img1, img0, img0));
figure; imagesc(rgb); colormap(gray(256)); axis image; colorbar;

%figure; imagesc(rgbOverlay(img, matchedMask, [1 0 0])); colormap(gray(256)); axis image; colorbar;

cellMask = [];



function newCC = updateEdgeInfo(matchedMask, CC)

labels = double(labelmatrix(CC));

newCC = bwconncomp(matchedMask, 8);
newCC.NumPixels = cellfun(@numel, newCC.PixelIdxList);
labelMap = mat2cell(labels(vertcat(newCC.PixelIdxList{:})), newCC.NumPixels,1);
labelMap = cellfun(@unique, labelMap, 'UniformOutput', false);

N = newCC.NumObjects;
% newCC.PixelOrder = cell(N,1);
newCC.EndpointIdx = cell(N,1);
newCC.rval = cell(N,1);
newCC.lval = cell(N,1);
for k = 1:N
    idx = setdiff(labelMap{k}, 0);
    %newCC.PixelOrder{k} = vertcat(CC.PixelOrder{idx});
    newCC.EndpointIdx{k} = horzcat(CC.EndpointIdx{idx});
    newCC.rval{k} = vertcat(CC.rval{idx});
    newCC.lval{k} = vertcat(CC.lval{idx});
end





% %------------------------------------------------------------------------------
% % III. Join remaining segments/endpoints using graph matching
% %------------------------------------------------------------------------------
% 
% % Remove long spurs
% cellBoundary = bwmorph(cellBoundary, 'thin');
% cellBoundary = bwmorph(cellBoundary, 'spur', 100);
% cellBoundary = bwmorph(cellBoundary, 'clean'); % spur leaves single pixels -> remove
% 
% % Create mask, use largest connected component within coarse threshold (removes potential loops in boundary)
% maskCC = bwconncomp(~cellBoundary, 4);
% csize = cellfun(@(c) numel(c), maskCC.PixelIdxList);
% [~,idx] = sort(csize, 'descend');
% % two largest components: cell & background
% int1 = mean(img(maskCC.PixelIdxList{idx(1)}));
% int2 = mean(img(maskCC.PixelIdxList{idx(2)}));
% cellMask = zeros(ny,nx);
% if int1 > int2
%     cellMask(maskCC.PixelIdxList{idx(1)}) = 1;
% else
%     cellMask(maskCC.PixelIdxList{idx(2)}) = 1;
% end
% 
% % loop through remaining components, check whether part of foreground or background
% for i = idx(3:end)
%     px = coarseMask(maskCC.PixelIdxList{i});
%     if sum(px) > 0.6*numel(px)
%         cellMask(maskCC.PixelIdxList{i}) = 1;
%     end
% end
% cellMask = imdilate(cellMask, strel('disk',1));
% 
% % Optional: erode filopodia-like structures
% if ~isempty(ip.Results.RemoveRadius)
%     cellMask = imopen(cellMask, strel('disk', ip.Results.RemoveRadius));
% end
%     
% % Final contour: pixels adjacent to mask
% B = bwboundaries(cellMask);
% cellBoundary = zeros(ny,nx);
% cellBoundary(sub2ind([ny nx], B{1}(:,1), B{1}(:,2))) = 1;
% 
% 
% 
% function out = connectEndpoints(inputPoints, queryPoints, radius, labels, cellBoundary, updateBoundary)
% if nargin<6
%     updateBoundary = true;
% end
% 
% dims = size(cellBoundary);
% out = zeros(dims);
% nq = size(queryPoints,1);
% [idx, dist] = KDTreeBallQuery(inputPoints, queryPoints, radius);
% 
% labSelf = labels(sub2ind(dims, queryPoints(:,2), queryPoints(:,1)));
% labAssoc = cellfun(@(i) labels(sub2ind(dims, inputPoints(i,2), inputPoints(i,1))), idx, 'UniformOutput', false);
% 
% % idx of endpoints belonging to other edges
% otherIdx = arrayfun(@(i) labAssoc{i}~=labSelf(i), 1:nq, 'UniformOutput', false);
% 
% % remove segment self-association (and thus query self-association)
% idx = arrayfun(@(i) idx{i}(otherIdx{i}), 1:nq, 'UniformOutput', false);
% dist = arrayfun(@(i) dist{i}(otherIdx{i}), 1:nq, 'UniformOutput', false);
% 
% % generate edge map
% E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:nq, 'UniformOutput', false);
% E = vertcat(E{:});
% 
% if ~isempty(E)
%     idx = E(:,1) < E(:,2);
%     
%     E = E(idx,:); % remove redundancy
%     
%     % generate weights
%     D = vertcat(dist{:});
%     D = D(idx);
%     
%     D = max(D)-D;
%     M = maxWeightedMatching(size(inputPoints,1), E, D);
%     
%     E = E(M,:);
%     
%     % add linear segments corresponding to linked endpoints
%     for i = 1:size(E,1)
%         iseg = bresenham([queryPoints(E(i,1),1) queryPoints(E(i,1),2)],...
%             [inputPoints(E(i,2),1) inputPoints(E(i,2),2)]);
%         out(sub2ind(dims, iseg(:,2), iseg(:,1))) = 1;
%     end
% end
% 
% if updateBoundary
%     out = double(out | cellBoundary);
% end
