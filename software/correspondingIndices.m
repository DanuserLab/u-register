function iCorrespond = correspondingIndices(pointsIn,curve2,maxDist)

if nargin < 3 || isempty(maxDist)
    maxDist = Inf;
end

%Finds the indices for the points in curve 2 which are closest to the input
%points
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

nPoints = size(pointsIn,2);


for j = 1:nPoints
   
    %Calculate distance to each point on curve 2
    currDists = sqrt( (curve2(1,:) - pointsIn(1,j)).^2 + (curve2(2,:) - pointsIn(2,j)).^2 );
    
    %Find closes point
    [minDist,iTmp] = min( currDists );
    if minDist <= maxDist
        iCorrespond(j) = iTmp;
    else
        iCorrespond(j) = NaN;
    end
    
end