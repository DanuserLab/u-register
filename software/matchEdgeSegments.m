% Francois Aguet, 10/2012
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

function [matchLabel, matchIndex] = matchEdgeSegments(CC, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('SearchRadius', 10, @isscalar);
ip.addParamValue('Mode', 'KSDistance', @(x) any(strcmpi(x, {'KSDistance', 'Edge'})));
ip.parse(varargin{:});

dims = CC.ImageSize;
labels = double(labelmatrix(CC));

pidx = vertcat(CC.PixelIdxList{:});

endpointIdx = vertcat(CC.EndpointIdx{:});
endpointLabels = labels(endpointIdx);

[yi, xi] = ind2sub(dims, pidx);
X = [xi yi];
[yi, xi] = ind2sub(dims, endpointIdx);
% use endpoints to find closest points on other edges
[idx, dist] = KDTreeBallQuery(X, [xi yi], ip.Results.SearchRadius);

np = numel(endpointIdx);
matchLabel = NaN(np,2);
matchIndex = NaN(np,2);
matchDist = NaN(np,1);
% loop through endpoints, determine matches
for k = 1:np
    % remove self-queries
    rmIdx = labels(pidx(idx{k}))==endpointLabels(k);
    idx{k}(rmIdx) = [];
    dist{k}(rmIdx) = [];
    
    % label of closest edge %%%% USE STRONGEST EDGE INSTEAD!!!
    % pick edge with highest intensity among closest 3
    if ~isempty(idx{k})
        iAvgInt = CC.AvgInt(labels(pidx(idx{k})));
        
        %if iAvgInt==max(iAvgInt)
        % match: first 'idx' with max 'iAvgInt'
        li = min(3, numel(iAvgInt));
        mi = find(iAvgInt(1:li)==max(iAvgInt(1:li)), 1, 'first');
        %mi = 1;
%         if mi<4
            % label of the matched points
            iMatchLabel = [endpointLabels(k) labels(pidx(idx{k}(mi)))];
            [~,sidx] = sort(iMatchLabel);
            matchLabel(k,:) = iMatchLabel(sidx);
            % index of the matched points
            iMatchIndex = [endpointIdx(k) pidx(idx{k}(mi))];
            matchIndex(k,:) = iMatchIndex(sidx);
            % distance btw matched points
            matchDist(k) = dist{k}(mi);
%         end
    end
end
% determine unique pairs of edges (shortest distance)
[~,sidx] = sort(matchDist);
matchLabel = matchLabel(sidx,:);
matchIndex = matchIndex(sidx,:);
rmIdx = isnan(matchLabel(:,1));
matchLabel(rmIdx,:) = [];
matchIndex(rmIdx,:) = [];
[~,sidx] = sort(matchLabel(:,1));
matchLabel = matchLabel(sidx,:);
matchIndex = matchIndex(sidx,:);

N = size(matchLabel,1);
% POSSIBLY LEAVE OUT: 
valid = false(N,1);
valid(1) = true;
for k = 2:N
    valid(k) = ~any(matchLabel(1:k-1,1)==matchLabel(k,1) & matchLabel(1:k-1,2)==matchLabel(k,2));
end
matchLabel = matchLabel(valid,:);
matchIndex = matchIndex(valid,:);
N = size(matchLabel,1);

%%
tmp = double(labels~=0);
tmp(endpointIdx(:)) = 2;
figure; imagesc(tmp); colormap(gray(256)); axis image; colorbar;
hold on;
[yi1, xi1] = ind2sub(CC.ImageSize, matchIndex(:,1));
[yi2, xi2] = ind2sub(CC.ImageSize, matchIndex(:,2));
X = [xi1 xi2];
Y = [yi1 yi2];
plot(X',Y', 'b')
% return
%%

kval = zeros(CC.NumObjects,1);
for k = 1:CC.NumObjects
    [~,~,kval(k)] = kstest2(CC.nvalDist{k}(:), CC.pvalDist{k}(:));
end

switch ip.Results.Mode
    case 'KSDistance'
        % cost based on KS distance
        cost = zeros(N,1);
        for k = 1:N
            [~,~,ksLL] = kstest2(CC.nvalProx{matchLabel(k,1)}(:), CC.nvalProx{matchLabel(k,2)}(:));
            [~,~,ksHH] = kstest2(CC.pvalProx{matchLabel(k,1)}(:), CC.pvalProx{matchLabel(k,2)}(:));
            
            %[~,~,ksLL] = kstest2(CC.nvalDist{matchLabel(k,1)}(:), CC.nvalDist{matchLabel(k,2)}(:));
            %[~,~,ksHH] = kstest2(CC.pvalDist{matchLabel(k,1)}(:), CC.pvalDist{matchLabel(k,2)}(:));

            % penalize cost when left/right distributions of a segment are close
            %cost(k) = (1-max([ksLL ksHH])) * kval(matchLabel(k,1)) * kval(matchLabel(k,2)) * CC.AvgInt(matchLabel(k,1)) * CC.AvgInt(matchLabel(k,2));
            
            % similar intensities, penalize difference in distributions
            %cost(k) = 1-abs(CC.AvgInt(matchLabel(k,1))-CC.AvgInt(matchLabel(k,2)))...
            
            
            % reward linking of high intensity segments, penalize difference in distributions
            cost(k) = min(CC.AvgInt(matchLabel(k,1)),CC.AvgInt(matchLabel(k,2)))...
                   .* (1-max([ksLL ksHH]));
                   %- max([ksLL ksHH]);% * kval(matchLabel(k,1)) * kval(matchLabel(k,2));
            
            
            %cost(k) = min(CC.AvgInt(matchLabel(k,1)),CC.AvgInt(matchLabel(k,2)))...
            %       - max([ksLL ksHH]);% * kval(matchLabel(k,1)) * kval(matchLabel(k,2));
            
            %cost(k) = cost(k) * min(ks1,ks2);
            %cost(k) = cost(k) * (1-(1-ks1)*(1-ks2));
            
        end
%         rmIdx = cost<0.2;
%         matchLabel(rmIdx,:) = [];
%         matchIndex(rmIdx,:) = [];
%         cost(rmIdx) = [];
    case 'Edge'
        % score edges: background vs. foreground
        medDiff = cellfun(@(i) nanmedian(i(:)), CC.rval) - cellfun(@(i) nanmedian(i(:)), CC.lval);
        medDiff = (medDiff-min(medDiff)) / (max(medDiff)-min(medDiff));
        nsize = CC.NumPixels/sum(CC.NumPixels); % ! geodesic distance should be used here !
        % cost = mean(medDiff(matchList),2)
        % reward strong and long edges
        cost = sum(nsize(matchLabel),2).*mean(medDiff(matchLabel),2);
        %rmIdx = cost<0.1;
        % matchList(rmIdx,:) = [];
        % cost(rmIdx) = [];
end
% figure; plot(cost);
% return

% [yi1, xi1] = ind2sub(CC.ImageSize, matchIndex(:,1));
% [yi2, xi2] = ind2sub(CC.ImageSize, matchIndex(:,2));
% X = [xi1 xi2];
% Y = [yi1 yi2];
% plot(X',Y', 'r--')
% for k = 1:numel(cost)
%    text(mean(xi1(k),xi2(k)), mean(yi1(k),yi2(k)), num2str(cost(k), '%.2f'), 'Color', 'r');
% end


M = maxWeightedMatching(CC.NumObjects, matchLabel, cost); % returns index (M==true) of matches
matchLabel = matchLabel(M,:);
matchIndex = matchIndex(M,:);

[yi1, xi1] = ind2sub(CC.ImageSize, matchIndex(:,1));
[yi2, xi2] = ind2sub(CC.ImageSize, matchIndex(:,2));
X = [xi1 xi2];
Y = [yi1 yi2];
plot(X',Y', 'r--')
% cost = cost(M);
% for k = 1:numel(cost)
%    text(mean(xi1(k),xi2(k)), mean(yi1(k),yi2(k)), num2str(cost(k), '%.2f'), 'Color', 'r');
% end



