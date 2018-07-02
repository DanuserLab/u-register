%nms = nonMaximumSuppression(res, th, sup)
%
% Inputs:   res : response
%            th : orientation
%           sup : value to set suppressed pixels to
%
% Uses the grid conventions of steerableDetector()
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

% Francois Aguet

function res = nonMaximumSuppression(res, th, sup)

if(nargin < 3)
    sup = 0;
end

[ny,nx] = size(res);

res = padarrayXT(res, [1 1], 'symmetric');

[x,y] = meshgrid(1:nx,1:ny);

% +1 interp
A1 = interp2(res, x+1+cos(th), y+1+sin(th),'linear',0);

% -1 interp
A2 = interp2(res, x+1-cos(th), y+1-sin(th),'linear',0);

res = res(2:end-1,2:end-1);

res(res<A1 | res<A2) = sup;
