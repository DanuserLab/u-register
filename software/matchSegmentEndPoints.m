%[mask, matchMask] = matchSegmentEndPoints(mask, theta, varargin) refines a binary edge/ridge map using graph matching to close gaps
%
% Inputs: 
% 
%    mask : binary edge/ridge map (output from steerable filter or similar)
%   theta : orientation map (from steerable filter or similar)
%
% Output:
%
%   matchedMask : refined mask
%  unmatchedIdx : index of unmatched segment endpoints
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

% Francois Aguet, 01/22/2012 (last modified 10/25/2012)

function [matchedMask, unmatchedIdx] = matchSegmentEndPoints(CC, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('CC', @isstruct);
ip.addParamValue('SearchRadius', 4, @isscalar);
ip.addParamValue('KeepJunctions', true, @islogical);
ip.addParamValue('Display', false, @islogical);
ip.parse(CC, varargin{:});
R = ip.Results.SearchRadius;

%=============================================
% I. Match segments
%=============================================

iter = 1;
matchesFound = true;
while matchesFound
    
    % index of all endpoints
    endpointIdx = [CC.EndpointIdx{:}]';
    % corresponding edge label
    endpointLabel = repmat(1:CC.NumObjects, [2 1]);
    endpointLabel = vertcat(endpointLabel(:));
    
    theta = cellfun(@(i) i([1 end]), CC.Theta, 'unif', 0);
    theta = vertcat(theta{:});
    
    [ye, xe] = ind2sub(CC.ImageSize, endpointIdx);    
    X = [xe ye];
    
    [idx, dist] = KDTreeBallQuery(X, X, R);
    
    if ~isempty(idx)
        
        % Generate all pairs that resulted from the query.
        % These pairs are the edges in the graph used for matching
        E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i} dist{i}], 1:numel(endpointIdx), 'UniformOutput', false);
        E = vertcat(E{:});
        
        % remove redundant pairs & self-queried points
        E = E(E(:,1) < E(:,2),:);
        
        % sort endpoints by edge label and distance
        % if multiple matches between two segments, retain connection with shortest distance
        M = sortrows([endpointLabel(E(:,1:2)) E(:,3) E(:,1:2)]);
        [~,ia] = unique(M(:,1:2), 'rows', 'first');
        M = M(ia,:);
        E = M(:,4:5);
        
        % Omitting this in the first pass adds stability by essentially ignoring (self matching)
        % small segments adjacent to (and offset from) gaps
        if iter>1
            % remove queries on same segment
            E(endpointLabel(E(:,1))==endpointLabel(E(:,2)),:) = [];
        end
        
        % consider only matches to the X highest-intensity segments
        T = sort(CC.AvgRes);
        T = T(end-10);
        E(CC.AvgRes(endpointLabel(E(:,1)))<T & CC.AvgRes(endpointLabel(E(:,2)))<T,:) = [];
        
        % Matching criteria:
        % 1) angle between the two endpoint vectors (orientations)
        t1 = theta(E(:,1));
        t2 = theta(E(:,2));
        a1 = abs(t1 - t2);
        a1(a1>pi) = a1(a1>pi)-pi;
        a2 = abs(a1-pi);
        minAngle = min(a1,a2);
        % cost function in [-1,1]: consider angles up to pi/4, penalize larger values
        cost = cos(2*minAngle);
        
        % 2) angle between the vector connecting the endpoints and the average endpoint orientation
        % endpoint vectors
        v1 = [cos(t1) sin(t1)]';
        v2 = [cos(t2) sin(t2)]';
        % angle btw. vectors
        dt = acos(sum(v1.*v2,1));
        v2(:,dt>pi/2) = -v2(:,dt>pi/2); % flip if wrong direction
        % vectors between endpoint pairs
        vL = [X(E(:,2),1)-X(E(:,1),1) X(E(:,2),2)-X(E(:,1),2)]';
        vL = vL./repmat(sqrt(sum(vL.^2,1)), [2 1]);
        % mean btw. endpoint vectors
        vMean = (v1+v2) ./ repmat(sqrt(sum((v1+v2).^2,1)), [2 1]);
        diffT = acos(sum(vMean.*vL,1))';
        diffT(diffT>pi/2) = pi-diffT(diffT>pi/2);
        % cost function: differences up to pi/4 are allowed
        cost(diffT>=pi/4) = -1;
        cost(diffT<pi/4) = cost(diffT<pi/4) .* (1-2*sin(2*diffT(diffT<pi/4)).^8);
        
        M = maxWeightedMatching(numel(endpointIdx), E, cost); % returns index (M==true) of matches
%%
        if ip.Results.Display && any(M)
            tmp = zeros(CC.ImageSize);
            for k = 1:CC.NumObjects
                tmp(CC.PixelIdxList{k}) = CC.AvgRes(k);
            end
            
            figure; imagesc(tmp); colormap(gray(256)); axis image; colorbar;
            hold on;
            % all endpoints
            %plot(xe(unmatchedIdx), ye(unmatchedIdx), 'rx');
            %T = theta(endpointIdx);
            % endpoint orientations
            quiver(X(:,1), X(:,2), cos(theta), sin(theta),0);            
            % endpoint candidates for matching
            %plot(X(unique(E(:)),1), X(unique(E(:)),2), 'go');
            
            plot([X(E(:,1),1) X(E(:,2),1)]', [X(E(:,1),2) X(E(:,2),2)]', 'y');
            mux = (X(E(:,1),1)+X(E(:,2),1))/2;
            muy = (X(E(:,1),2)+X(E(:,2),2))/2;
            %for k = 1:numel(mux)
            %   text(mux(k)+0.2, muy(k), num2str(cost(k), '%.2f'), 'Color', 'c', 'VerticalAlignment', 'bottom')
            %end
            
            % plot mean vector for each pair
            quiver(mux, muy, vMean(1,:)', vMean(2,:)',0, 'g');
            
            % plot transverse vector for each pair
            %quiver(mux, muy, vL(1,:)', vL(2,:)',0, 'r');
            %quiver(mux, muy, -vL(1,:)', -vL(2,:)',0, 'r');
            %quiver(mux, muy, vL(2,:)', -vL(1,:)',0, 'r');
            %quiver(mux, muy, vL(2,:)'-vL(1,:)', -vL(1,:)'-vL(2,:)',0, 'r');

            %plot([X(E(rmIdx,1),1) X(E(rmIdx,2),1)]', [X(E(rmIdx,1),2) X(E(rmIdx,2),2)]', 'g--');
            plot([X(E(M,1),1) X(E(M,2),1)]', [X(E(M,1),2) X(E(M,2),2)]', 'r');
            title(['Iteration ' num2str(iter)]);
        end
        %%
        % retain only pairs that are matches
        E = E(M,:);
        
        % remove 'self' matches (only relevant for 1st iteration)
        if iter==1
%             E(endpointLabel(E(:,1))==endpointLabel(E(:,2)),:)
            E(labels(endpointIdx(unmatchedIdx(E(:,1))))==labels(endpointIdx(unmatchedIdx(E(:,2)))), :) = [];
        end
        
        % update CC with matches
        
        % fill mask
        for i = 1:size(E,1)
            iseg = bresenham([X(E(i,1),1) X(E(i,1),2)], [X(E(i,2),1) X(E(i,2),2)]);
            matchedMask(sub2ind([ny nx], iseg(:,2), iseg(:,1))) = 1;
        end
        
        matchedIdx = unmatchedIdx(unique(E(:))); % current iter only
        unmatchedIdx = setdiff(unmatchedIdx, matchedIdx);
    end
    
    iter = iter + 1;
    if size(E,1)==0
        matchesFound = false;
    end
end

matchedMask = segmentMatrix | matchedMask;
if ip.Results.KeepJunctions
    matchedMask = matchedMask | junctionMatrix;
end

unmatchedIdx = endpointIdx(unmatchedIdx);
