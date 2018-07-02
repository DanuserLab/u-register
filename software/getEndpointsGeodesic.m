% Francois Aguet, 10/28/2012
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

function endpointIdx = getEndpointsGeodesic(segmentMask, endpointIdx)

segmentMask = double(segmentMask~=0);

% if no endpoints are provided, find pixels with a single neighbor
if nargin<2 || isempty(endpointIdx)
    nn = double(bwmorph(segmentMask, 'thin'));
    nn = (imfilter(nn, ones(3), 'same')-1) .* nn;
    endpointIdx = find(nn==1);
end
endpointIdx = endpointIdx(:);

CC = bwconncomp(segmentMask, 8);
labels = double(labelmatrix(CC));

% assign endpoints to connected components, and sort lists
endpointLabels = labels(endpointIdx);
[endpointLabels, sidx] = sort(endpointLabels);
tmp = diff([0 find([diff(endpointLabels(:)')~=0 1])])';
nEndpoints = zeros(CC.NumObjects,1);
nEndpoints(unique(endpointLabels)) = tmp; % in case there are segments w/o endpoints
endpointIdx = endpointIdx(sidx);


% initialize array to store endpoint-to-furthest endpoint distances
endpointDist = zeros(size(endpointIdx));
segmentMask = logical(segmentMask);
cumidx = [0; cumsum(nEndpoints(1:end-1))];
for n = 1:max(nEndpoints)
    % nth endpoint for each CC
    nidx = cumidx(nEndpoints>=n)+n;
    D = bwdistgeodesic(segmentMask, endpointIdx(nidx));
    
    % store max. distance
    endpointDist(nidx) = cellfun(@(i) max(D(i)), CC.PixelIdxList(nEndpoints>=n));
end

% retrieve 2 endpoints with largest distances for each CC
endpointDist = mat2cell(endpointDist, nEndpoints, 1);
endpointIdx = mat2cell(endpointIdx, nEndpoints, 1);
for c = 1:CC.NumObjects
    [~,idx] = sort(endpointDist{c}, 'descend');
    endpointIdx{c} = endpointIdx{c}(idx(1:min(numel(idx),2)));
end
endpointIdx = endpointIdx';
% CC.EndpointIdx = endpointIdx';
