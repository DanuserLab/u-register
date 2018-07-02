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

function [edgeMask] = linkEdgeSegments(edgeMask, linkIndex)

edgeMask = edgeMask~=0;
dims = size(edgeMask);
linkMask = zeros(dims);
for k = 1:size(linkIndex,1)
    [y0, x0] = ind2sub(dims, linkIndex(k,1));
    [y1, x1] = ind2sub(dims, linkIndex(k,2));
    iseg = bresenham([x0 y0], [x1 y1]);
    iseg = sub2ind(dims, iseg(:,2), iseg(:,1));
    linkMask(iseg) = 1;
end
linkMask = bwmorph(linkMask, 'thin');
edgeMask = edgeMask | linkMask;
