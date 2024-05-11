function closedContours = closeContours(contoursIn,maskIn,shiftVal)
%CLOSECONTOURS takes the input  open contours and closes them by filling in the gaps where they meet the image border
%
% closedContours = closeContours(contoursIn,maskIn)
% closedContours = closeContours(contoursIn,maskIn,shiftVal)
%
% Description:
% 
% The contours as returned by the matlab contouring functions (contour.m,
% contourc.m etc) will contain gaps where a given contour meets the edge of
% the matrix is being contoured. This function closes these contours by
% filling in the appropriate areas at the image boundary, so that the
% closed contour encloses only areas of the mask with true (nonzero) values.
% 
% Required Input: 
% 
%   contoursIn - A 1xM cell array of the M input contours, as returned by
%                separateContours.m
% 
%   maskIn -    The mask the contours were derived from. 
% 
% Optional Input:
% 
%   shiftVal - If greater than zero, the closed contour values at the
%   matrix border will be shifted this far outside the matrix. Optional.
%   Default is 0 (no shift). Must be less than .5.
%
% Output:
% 
%   closedContours - A 1xM array of the closed contours. 
% 
% 
% Hunter Elliott
% 4/2010
% 
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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


%% ----- Input ----- %%

if nargin < 2 || isempty(maskIn) || isempty(contoursIn)
    error('Must input both a cell-array of contours and the matrix which the contours came from!')
end

if ~iscell(contoursIn)
    error('The input contours must be separated into a cell array! Try using separateContours.m!')
end

if nargin < 3 || isempty(shiftVal)
    shiftVal = 0;
elseif abs(shiftVal) >= .5
    error('Shiftval must be smaller than .5!')
end


%% -----  Init ----- %%

[M,N] = size(maskIn);

nContours = length(contoursIn);

%Get the coordinates of the image border in a  clock-wise fashion, starting
%at origin.
if shiftVal ==0%Do this as an if/else to help readability    
    borderCoord = vertcat([1:(N-1)      ones(1,M-1)*N  N:-1:2         ones(1,M-1)],...
                          [ones(1,N-1)  1:(M-1)        ones(1,N-1)*M  M:-1:2]);
else
    %Do this as an if/else to help readability
    borderCoord = vertcat([1:N              ones(1,M)*N+shiftVal  N:-1:1                  ones(1,M)-shiftVal],...
                          [ones(1,N)-shiftVal  1:(M)                ones(1,N)*M+shiftVal  M:-1:1]);
end
    
closedContours = cell(nContours,1);

%% ----- Closure ---- %%

for j = 1:nContours
    
    %Verify that this contour needs closure by checking that the first
    %point touches the image border.    
    iTouchStart = find(all(bsxfun(@eq,round(borderCoord),round(contoursIn{j}(:,1)))),1,'first');    
    
    if ~isempty(iTouchStart)
        %Find where the last point touches the border        
        iTouchEnd = find(all(bsxfun(@eq,round(borderCoord),round(contoursIn{j}(:,end)))),1,'last');
        
        %Close the curve in both directions, then compare the values on the
        %two closed curves. 
        %Previous attempts to choose between the closure directions based
        %on the values only very near the start/end point failed due to the
        %discretized nature of the input matrix and some quirks of the
        %contours produced by contourc.m

        %Make sure the curve touches the border at more than one point.
        if iTouchEnd ~= iTouchStart
        
            %Check if the border used for closing crosses the origin
            if iTouchEnd < iTouchStart

                %Close clockwise, adding border points to the contour in their original order
                closeClockwise = borderCoord(:,iTouchEnd+1:iTouchStart-1);
                
                %Close counter-clockwise, adding border points to the
                %contour in reverse order and looping past origin
                closeCounterClockwise = [borderCoord(:,iTouchEnd-1:-1:1) borderCoord(:,end:-1:iTouchStart+1)];

            elseif iTouchEnd > iTouchStart

                %Close clockwise, adding border points to the contour in their
                %original order and looping past the origin
                closeClockwise = [borderCoord(:,iTouchEnd+1:end) borderCoord(:,1:iTouchStart-1)];

                %Close counter-clockwise, adding the border points in reverse
                %order
                closeCounterClockwise = borderCoord(:,iTouchEnd-1:-1:iTouchStart+1);
            end
        
            %Sample the input matrix at the border coordinates for the two
            %closure directions to check which only contains true mask
            %areas.            
            
            if nnz(~maskIn(sub2ind([M,N],round(closeClockwise(2,:)),round(closeClockwise(1,:)))))==0;
                closedContours{j} = [contoursIn{j} closeClockwise];
            elseif nnz(~maskIn(sub2ind([M,N],round(closeCounterClockwise(2,:)),round(closeCounterClockwise(1,:)))))==0;
                closedContours{j} = [contoursIn{j} closeCounterClockwise];
            else
                error(['Check input matrix and contours - unable to close contour ' num2str(j) ])
            end                                        
            
        else
            
            %No need to do anything - this curve just touches the border
            %at one point.
            closedContours{j} = contoursIn{j};            
            
        end
    else    
        %No need to do anything - this curve doesn't touch the border
        closedContours{j} = contoursIn{j};
    end
        
end

