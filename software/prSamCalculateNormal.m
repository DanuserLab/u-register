function nxy_n=prSamCalculateNormal(img, sp_n, r_n)
% prSamCalculateNormal calculates the unit normals relative to the edge 
%               This function calculates the unit normals along a edge
%               given by the two component spline at the paremeters given
%               by "r_n". The normal points away from the object assuming
%               that the object has a higher average intensity.
%
% SYNOPSIS    nxy_n=prGetProtRegion(img, sp_n, r_n)
%
% INPUT       img       : normalized grayscale image
%             sp_n      : B-spline 
%             r_n       : spline parameter
%
% OUTPUT      nxy_n      : edge unit normal vector
%   
% First Created by Matthias Machacek 11/11/03 (prGetProtRegion.m)
% Modefied as prSamCalculateNormal.m by Shann-Ching Chen, LCCB, 08/26/2008
% See also: prSamProtrusion
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


[n_img, m_img]=size(img);
knots_nr=length(r_n);


xy_n=fnval(sp_n,r_n)';
x_n=xy_n(:,1);  y_n=xy_n(:,2);

%derivative spline 
sp_dxy_n = fnder(sp_n);

%derivatives at discrite locations
dxy_n = fnval(sp_dxy_n,r_n)';
dx_n=dxy_n(:,1);
dy_n=dxy_n(:,2);
%normalize
l=sqrt(dx_n.^2+dy_n.^2);
dx_nn=dx_n./l;
dy_nn=dy_n./l;
%the normal unit vector (not oriented!!)
nx_n= dy_n;
ny_n=-dx_n;

%determine the object side of the edge
p1_x=round(x_n+2*nx_n);
p1_y=round(y_n+2*ny_n);

p2_x=round(x_n-2*nx_n);
p2_y=round(y_n-2*ny_n);

in_out1=0;
in_out2=0;

for i=1:length(p1_x)
    logic_exp1 = p1_x(i) >= 1     & p1_y(i) >= 1     & p2_x(i) >= 1     & p2_y(i) >= 1;
    logic_exp2 = p1_x(i) <= m_img & p1_y(i) <= n_img & p2_x(i) <= m_img & p2_y(i) <= n_img;
    if logic_exp1 & logic_exp2
        in_out1 = in_out1 + double(img(p1_y(i), p1_x(i)));
        in_out2 = in_out2 + double(img(p2_y(i), p2_x(i)));
    end
end

%assume that the background has a lower intensity
%the normal is pointin away from the object!
if in_out1 > in_out2
   nx_n= -nx_n;
   ny_n= -ny_n; 
end

nxy_n = [nx_n ny_n];
