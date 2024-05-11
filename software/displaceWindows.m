function dispWindows = displaceWindows(windowsIn,vX,vY,edget1,edget2,protrusion)
%DISPLACEWINDOWS moves the input windows based on the input vector field
%
% dispWin = displaceWindows(windowsIn,X,Y,edge,protrusion);
%
% Displaces the input windows based on the input vector field and
% protrusion vectors. The output windows are moved by an amount equal to
% the value of the vector field interpolated to the position of the input
% window vertices.
%
%
% Input:
% 
%   windowsIn - The cell array containing the windows to displace, as
%   created using getMaskWindows.m
%
%   [X,Y] - 2D matrices specifying the X and Y components respectively of
%   the vector field at each pixel in the windowed image.
%
%   edget1 - 2xM vector containing the XY positions of the object border in
%   the current frame
%
%   edget2 - 2xN vector containing the XY positions of the object border in
%   the next frame
%
%   protrusion - A 2xM vector with the protrusion vectors at the border of
%   the object which has been windowed.
% 
% Output:
% 
%   dispWin - The displaced windows.
%
%
% Hunter Elliott
% Re-written 8/2010
%
%% ------ Input ------- %%
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


if nargin < 6 || isempty(windowsIn) || isempty(vX) || isempty(vY) ...
        || isempty(edget1) || isempty(edget2) || isempty(protrusion)
    error('You must input all 5 variables!');
end


%% ------ Init ----- %%

nStrips = numel(windowsIn);
dispWindows = cell(1,nStrips);

[mVec,nVec] = size(vX);
[mT,nT] = size(vY);

if mT ~= mVec || nT ~= nVec
    error('X & Y vector fields must be same size!');        
end

[X,Y] = meshgrid(1:nVec,1:mVec);
   


%% --------- Displace --------- %%
    
%Loop through windows and displace.
for iStrip = 1:nStrips
    
    nBands = numel(windowsIn{iStrip});
    dispWindows{iStrip} = cell(1,nBands);
    
    for iBand = 1:nBands

        %NOTE - THIS IS DUMB - You are repeating points where windows share edges!!!Save time by just doing outer borders (except on innermost
        %band)??? - HLE                

        %Faster to only send local area of window? Or to use
        %TriScatteredInterp and evaluate function multiple times? -HLE

        %Get interpolated X & Y vector values for window borders.
        %Windows are in matrix coord so X and Y are reversed.            

        %We treat the window sides on the edge of the object
        %differently, so don't displace these using vector field.
        if ~isempty(windowsIn{iStrip}{iBand})
            if iBand == 1;               
                nSides = numel(windowsIn{iStrip}{iBand})-1;

                %Displace the outer points using the protrusion vectors
                %themselves. This ensures that they will remain on the
                %object periphery.

                %Find points on smoothed edge which correspond to these
                iProt = correspondingIndices(windowsIn{iStrip}{iBand}{end},edget1);

                %Displace using these vectors, and find the closest
                %points on the edge in the next frame.
                dispEdge = windowsIn{iStrip}{iBand}{nSides+1} + protrusion(:,iProt);                
                iNextFrame = correspondingIndices(dispEdge,edget2);
                dispWindows{iStrip}{iBand}{nSides+1} = edget2(:,iNextFrame);


            else
                nSides = numel(windowsIn{iStrip}{iBand});
            end

            for iSide = 1:nSides                    

                nPts = size(windowsIn{iStrip}{iBand}{iSide},2);

                %Initialize displacement vectors
                dispVecs = nan(2,nPts);                            

                dispVecs(1,:) = interp2(X,Y,vX,...
                    windowsIn{iStrip}{iBand}{iSide}(1,:),...
                    windowsIn{iStrip}{iBand}{iSide}(2,:),'*linear');
                dispVecs(2,:) = interp2(X,Y,vY,...
                    windowsIn{iStrip}{iBand}{iSide}(1,:),...
                    windowsIn{iStrip}{iBand}{iSide}(2,:),'*linear');                  


                %Deal with NaNs - these are usually points which are close to
                %the object boundary, and get propagated outside the object
                %due to numerical error in solution.
                isNan = any(isnan(dispVecs),1);                    
                if any(isNan)

                    %Better way to do this??? - HLE
                    windowsIn{iStrip}{iBand}{iSide}(:,isNan) = [];
                    dispVecs(:,isNan) = [];


                end

                %Displace window borders
                dispWindows{iStrip}{iBand}{iSide} = windowsIn{iStrip}{iBand}{iSide} + dispVecs;


            end
        end
    end
end

