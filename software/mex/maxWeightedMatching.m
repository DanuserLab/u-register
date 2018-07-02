function M = maxWeightedMatching(n,E,W)
% function M = maxWeightedMatching(n,E,W)
% maxWeightedMatching computes the maximum weighted matching of a graph 
% consisting of n nodes connected by the edges defined in E and weighted 
% with W. A matching is a subset of the edges defined in E. This matching 
% maximizes the sum of the weights of the edges in the subset so that each 
% node is connected to another node not more than once.
%
% Required Inputs:
% n              Scalar defining the number of nodes
% E              An m x 2 array representing m edges. Each line consists of 
%                two node indices. The order of the indices does not matter.
% W              An m x 1 array representing the weights of the m edges.
% 
% Optional Inputs:
%
% Outputs:
% M              A m x 1 logical array representing the matching. True
%                means the edge is part of the matching.
%
% Example: 
%                o - o
%                | X |
%                o - o
%
%                n = 4;
%                E = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
%                W = [2;6;7;8;5;3];
%                M = maxWeightedMatching(n,E,W)
%
%
% Sylvain Berlemont, Pascal Berard, October 2011
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
