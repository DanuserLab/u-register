function [contourArray,contourValues] = separateContours(contoursIn)
%EXTRACTCONTOURS seperates the individual contours in the countour matrix output from contour.m or contourc.m
% 
% contourArray = separateContours(contoursIn)
% [contourArray,contourValues] = separateContours(contoursIn)
% 
% Desctiption:
%
%   The matlab contouring functions (e.g. contour.m and countourc.m) return
%   all the contours as a single 2xM matrix. This function seperates these
%   contours into a cell-array, with each cell containing a single contour
%   as a 2xN matrix with the x and y cooridnates of each point on that
%   contour.
% 
% See Also: cleanUpContours.m 
%
% Input:
% 
%   contoursIn - A 2xM contour matrix, as returned by contourc.m, contour.m
%                etc, where M >= 2
% 
% 
% Output:
% 
%   contourArray - A 1xP cell-array, with each element containin a single
%                  contour. These contours are stored as a 2xN matrix, with
%                  the x-values in the top row (and the y values in the
%                  bottom row)
% 
% 
%   contourValues - The iso-values that correspond to each contour in the
%                   cell array. (Depending on the topology, there may be
%                   multiple contours at a single iso-value.)
%
% 
% 
% Hunter Elliott 
% 4/2010 
% 
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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

%% ------- Input ------- %%

if nargin < 1 || isempty(contoursIn) || size(contoursIn,1) ~=2 || ...
        size(contoursIn,2) < 2 || ndims(contoursIn) ~= 2
    error('The 1st and only input must be a 2xM contour matrix, where M>=2!! Check input!')
end


%% ----- Contour Separation ------ %%

currPos = 1;%Starting position in contour matrix.
iContour = 1; %Counter for number of contours

%Over-initialize the contour arrays, since at this point we don't know how
%many contours there will be.
nTotal = size(contoursIn,2); 
contourArray = cell(ceil(nTotal/2),1);
contourValues = zeros(ceil(nTotal/2),1);

while true

    if currPos < nTotal        
        %Get the size, in points, of the current contour. This is stored in the
        %contour matrix.
        currSize = contoursIn(2,currPos);
        
        %Use this size to store this contour in the array
        contourArray{iContour} = contoursIn(:,currPos+1:currPos+currSize);  
        
        %Get the contour value,
        contourValues(iContour) = contoursIn(1,currPos);
        
        %Advance our position in the matrix
        currPos = currPos + currSize + 1;
        
        %count this contour
        iContour = iContour+1;

    else
        %Unless we are at the last contour...
        break
    end

end

%% ------ Output ----- %%


%Remove over-initialized elements in arrays
contourArray = contourArray(1:iContour-1);
contourValues = contourValues(1:iContour-1);


