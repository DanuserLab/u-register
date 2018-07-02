function L = prSamProtrusionOptSp ( w )
% prSamProtrusionOpt is the funciton of new Protrusion model based
% on constrainted optimization
% 
% Input and Output: w is the solution for best spline parameter sets of pixels at time t
%
% See also: prSamProtrusion, prSamProtrusionOptSp
%
% Last updated: August 26, 2008 by Shann-Ching Chen, LCCB
% See also: prSamProtrusion, prSamProtrusionDownSpl
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

global D
%% optimize over w

L = 0;
spline_pixel=fnval(D.sp_s,w)';
D.iter = D.iter + 1;
D.L1(D.iter) = sum(sum((D.pixel_tm1 - spline_pixel).^2,2));
%D.L2(D.iter) = sum(sum(diff(spline_pixel).^2,2));
D.L2(D.iter) = sum(sum(diff(w).^2,2));

if D.iter == 1
    D.alpha = D.L1(D.iter)/D.L2(D.iter)/D.ratio;
    tmp = D.L1(D.iter) + D.alpha*D.L2(D.iter);
    %fprintf(1,'iter = %5d, L=%.2f(%.2f+%.2f)\n', D.iter, tmp, D.L1(D.iter),D.alpha*D.L2(D.iter));
    clear tmp;
end

L = D.L1(D.iter) + D.alpha*D.L2(D.iter);
%fprintf(1,'iter = %5d, L=%.2f(%.2f+%.2f)\n', D.iter, L, D.L1(D.iter),D.alpha*D.L2(D.iter));

