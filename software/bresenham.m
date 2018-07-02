function P = bresenham(p0, p1, conn)
%BRESENHAM computes the integer positions on a line between the
% positions xS and xE
%
% SYNOPSIS P=bresenham(p0,p1,conn)
%
% INPUT p0 : coordinate of line start point
%       p1 : coordinate of line end point
%     conn : connectivity of the line. The connectivity is either 4 or 8.
%            The default is conn=8 (if conn is not specified). For a
%            c4-connected line enter conn=4.
%
% OUTPUT P : nx2 matrix with the coordinates of all the
%            integer positions on the line
%
% Sylvain Berlemont, 2009
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

dx = p1(2) - p0(2);
dy = p1(1) - p0(1);

incr1 = [1 1];
incr2 = [0 0];

if dx < 0
    dx = -dx;
    incr1(2) = -1;
end

if dy < 0
    dy = -dy;
    incr1(1) = -1;
end

if dy >= dx
    dqr = 2 * dx;
    dqru = dqr - 2 * dy;
    q = dqr - dy;
    l = dy;
    incr2(1) = incr1(1);
else
    dqr = 2 * dy;
    dqru = dqr - 2 * dx;
    q = dqr - dx;
    l = dx;
    incr2(2) = incr1(2);
end

P = zeros(l + 1, 2);
p = p0;

for d = 1:l+1
    P(d, :) = p;
    
    if (q > 0)
        p = p + incr1;
        q = q + dqru;
    else
        p = p + incr2;
        q = q + dqr;
    end
end

% This is code added by Achim, 2010. Note that
% exchanging the order of p0 and p1 gives different results, because the
% P from above is not symmetric. The code that follows should give
% symmetric results:
if nargin>2 && conn==4
    n=size(P,1);
    i=1;
    while i<n
        % check if two consecutive points are connected through a vertex
        % instead of an edge:
        dxy=P(i+1,:)-P(i,:);
        if sum(abs(dxy))==2 % if there is a shift in x and y:
            % create new point with the x component of the first point and
            % the y component of the second point, or vice verca.
            % The if statement makes it symmetric (p0 <-> p1):
            if P(i,1)>P(i+1,1)
                xPt=[P(i,1)   P(i+1,2)];
            else
                xPt=[P(i+1,1)   P(i,2)];
            end
            
            % now insert this point:
            newP(1:i,:)=P(1:i,:);
            newP(i+1,:)=xPt;
            newP(i+2:n+1,:)=P(i+1:end,:);
            P=newP;
            clear newP;
            
            % the length has also changed now (important for the
            % termination of the loop):
            n=length(P);
        end
        i=i+1;
    end
end
