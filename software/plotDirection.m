function plotDirection(linesIn,varargin)
%PLOTDIRECTION plots the input 2D lines(s) with an arrow indicating the direction of the curve
%
% plotDirections(linesIn)
% plotContours(linesIn,plotStyle1,plotStyle2,...)
%
% Description:
%   
%   Plots the input 2D lines and shows the direction the points occur
%   within the line by placing a small arrow at the first point.  
% 
% Input:
% 
%   linesIn - A single 2xM or Mx2 matrix containing the line to plot, or a
%             cell array of 2xM or Mx2 matrices if multiple lines are to be
%             plotted.
% 
%   plotStyleString - Optional. A string (or strings) specifying the
%                     style/color to use when plotting the lines (same as
%                     with the plot command)
%
% Hunter Elliott
% 4/2010
% NOTE: Replaces my old plotContours.m function.
%
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

if ~iscell(linesIn)
    linesIn = {linesIn};
end

if nargin >  1
    plotArgs = varargin;
else
    plotArgs = {'b'};
end

hold on

%Convert all the lines to 2xM
needsTranspose = cellfun(@(x)(size(x,1) ~= 2),linesIn);
linesIn(needsTranspose) = cellfun(@(x)(x'),linesIn(needsTranspose),'UniformOutput',false);

%Plot the lines
cellfun(@(x)(plot(x(1,:),x(2,:),plotArgs{:})),linesIn)

%Plot the first point with arrow
iNotPt = cellfun(@(x)(size(x,2)> 1),linesIn); %Make sure the line has 2 pts!
cellfun(@(x)(quiver(x(1,1),x(2,1),x(1,2)-x(1,1),x(2,2)-x(2,1),10,plotArgs{:},'MaxHeadSize',10)),linesIn(iNotPt));
cellfun(@(x)(plot(x(1,1),x(2,1),plotArgs{:},'Marker','x')),linesIn(iNotPt));
