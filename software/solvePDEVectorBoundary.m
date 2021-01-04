function varargout = solvePDEVectorBoundary(xy,uv,pdePar,imgSize,meshQuality,nonLin)
%SOLVEPDEVECTORBOUNDARY solves the selected PDE using the input vectors as a boundary condition 
% 
%                     [X,Y] = solvePDEVectorBoundary(xy,uv)
%                 [u,p,e,t] = solvePDEVectorBoundary(xy,uv)
%                     ...   = solvePDEVectorBoundary(xy,uv,pdePar)
%                     ...   = solvePDEVectorBoundary(xy,uv,pdePar,imSize)
%                     ...   = solvePDEVectorBoundary(xy,uv,pdePar,imSize,meshQuality)
% 
% This function solves the selected partial differential equation (PDE)
% over the area enclosed by the vectors x and y given the vector-valued
% boundary conditions specified by the vectors u and v. The solution are
% found using the matlab PDE toolbox. The PDE must be of the form supported
% by the PDE toolbox function adaptmesh.m.
% 
% More specifically, the function finds the vector-valued solution u over
% the two dimensional domain omega defined by the input polygon xy, given
% the vector-valued Dirichlet boundary condition  u(omega) = uv. In the
% matlab PDE toolbox, vector valued solutions are treated as systems of
% equations. Therefore this corresponds to the solution of a system of
% elliptic PDEs
% 
% Input:
% 
%   xy - The 2xM or Mx2 matrix containing the x and y positions of the
%         boundary of the area to solve the PDE. These positions must form
%         a closed polygon (M>=3 and the last and first points are
%         adjacent).
% 
%   uv - The 2xM or Mx2 vectors containing the x and y components of the
%         boundary condition at each boundary point specified by xy.
%         Must be the same size as xy.
% 
%   pdePar - 1x3 Cell array specifying the parameters of the PDE to solve
%            in the object interior. The elements of this vector correspond
%            to the coefficients c, a and f used by the functions
%            assembpde.m and adaptmesh.m. See the help of assempde.m for
%            details on these coefficients and their specification.
%
%   imgSize - A 1x2 positive integer vector containing the size of the area
%             over which to solve the PDE. Areas outside of the area
%             specified by x and y will have zero values. Only used if X,Y
%             outputs are requested. Optional. If not specified, the X and
%             Y matrices will be just large enough to fit area enc)losed by
%             xy.
%
%   meshQualtiy - Integer scalar between 1-10. The quality of the
%                 triangular mesh to use when solving the PDE. Large
%                 numbers will DRASTICALLY increase computation time, while
%                 asymptotically decreasing error in the solution.
%                 Optional. Default is 3.
%
%   nonLin - 'off' or 'on' If 'on', the non-linear solver will be used.
%            This is required in cases where the PDE coefficients contain
%            the value of the solution OR it's first derivative. Optional.
%            Default is 'off' (use linear solver).
%
%
% Output:
%
%   [X,Y] - The rectangular matrices specifying the X and Y components of
%           the solution to the PDE at regularly spaced points within the
%           solution area. This output is slightly slower as it requires
%           that the solution defined on the triangular mesh be
%           interpolated to a homogeneous grid. 
%
%      ------------------------ OR -----------------------------
%
%   [u,p,e,t] -  The solution to the PDE on the triangular
%                mesh used by the PDE toolbox. This contains the solution
%                values uv, vertex locations p, edge topology e and
%                triangles t for the x and y components of the solution.
%                See assempde.m for details on their format.
%
% Hunter Elliott
% 8/2010
%
%
% Copyright (C) 2021, Danuser Lab - UTSouthwestern 
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

%% ----------- Input ---------- %%
if nargin < 2 || isempty(xy) || isempty(uv)
    error('Must input boundary location xy and boundary condition uv!')
end

if nargin < 3 || isempty(pdePar)
    error('You must specify the coefficients of the PDE you wish to solve!')
end

%Convert to Mx2 if input as 2xM
if size(xy,2) ~= 2
    xy = xy';
end
if size(uv,2) ~= 2
    uv = uv';
end

%Check that border is closed. 
if ~isequal(xy(1,:),xy(end,:))
    error('The input polygon defined by xy must be closed - that is, the first and last points must be identical!')
end

nPts = size(xy,1);

if nPts < 3 || size(uv,1) ~= nPts
    %NOTE: It is not actually required that xy and uv be the same length,
    %but it is an easy way to ensure correspondence between a location on
    %the boundary and the boundary condition.... HLE
    error('xy and uv must be matrices of the same size - Mx2 or 2xM where M >= 3!')
end

    
if nargin < 4
    imgSize = [];
end

if nargin < 5 || isempty(meshQuality)
    meshQuality = 3;
elseif (meshQuality <= 0) || (round(abs(meshQuality)) ~= meshQuality) || ...
    meshQuality > 10
    error('Mesh quality must be a positive integer between 1 and 10!')
end

if nargin < 6 || isempty(nonLin)
    nonLin = 'off';
end

%% --------- Init ----------- %%    

%We use global variables for the boundary coord and values, because the PDE
%toolbox is very strict about geometry and boundary condition inputs (it
%won't allow function handles or anonymous functions)
global BOUND_COND
global OBJ_BOUND

%Check if curve specified by xy runs clockwise or counterclockwise

%If it's not clockwise, we need to reverse it and the boundary conditions.
if ~isCurveClockwise(xy);
    xy = xy(:,end:-1:1);
    uv = uv(:,end:-1:1);
end

%Fit interpolating spline to boundary.
OBJ_BOUND = spline(linspace(0,1,nPts),xy');

OBJ_BOUND.ptsPerBS = 10;

%If the image size wasn't input, make sure the range includes the whole
%cell
if isempty(imgSize)
    sTmp = linspace(0,1,nPts);
    vs = ppval(OBJ_BOUND,sTmp);
    imgSize = ceil(max(vs,[],2) + 10); %Leave a little room...
end

%% -------- Solve ---------- %%
%Use adaptmesh to both refine the FEM mesh and to solve the PDE.

%Set global boundary condition variable to spline representation of uv
BOUND_COND = spline(linspace(0,1,nPts),uv');

%Initialize mesh
[p,e,t] = initmesh('boundaryGeometry');

Ui = 1;%Initial guess for solution.

%Do the specified number of rounds of adaptive refinement
[u,p,e,t] = adaptmeshHLE('boundaryGeometry','boundaryCondition',...
                         pdePar{1},pdePar{2},pdePar{3},...
                         'Ngen',meshQuality,'Mesh',p,e,t,...
                         'Nonlin',nonLin,'Init',Ui);
                     
                     
%% ------- Output ----- %%

if nargout > 2
    %If the solution was requested in mesh form
    varargout{1}=u;
    varargout{2}=p;
    varargout{3}=e;
    varargout{4}=t;    
else
    %Otherwise, interpolate the mesh solution to gridded data
    nT = size(p,2);
    %The PDE toolbox returns the x and y values in a single column vector
    %sequentially.
    varargout{1} = tri2grid(p,t,u(1:nT),1:imgSize(2),1:imgSize(1));
    varargout{2} = tri2grid(p,t,u(nT+1:end),1:imgSize(2),1:imgSize(1));
end



