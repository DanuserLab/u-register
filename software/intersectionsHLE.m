function [x0,y0,iout,jout] = intersectionsHLE(x1,y1,x2,y2,robust)
%INTERSECTIONSHLE Intersections of curves.
%
%   NOTE:Modified from version found on matlab exchange to fix a few bugs,
%   speed calculation, sort output etc. Part of this modification involves
%   removal of redundant intersections, and rounding of intersections
%   indices (iout, jout) slightly different from integer indices due
%   primarily to numerical error. If you need intersection indices which
%   are more precise than +/-1e-6, don't use this function. - HLE
%
%   Computes the (x,y) locations where two curves intersect.  The curves
%   can be broken with NaNs or have vertical segments.
%
% Example:
%   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% where X1 and Y1 are equal-length vectors of at least two points and
% represent curve 1.  Similarly, X2 and Y2 represent curve 2.
% X0 and Y0 are column vectors containing the points at which the two
% curves intersect.
%
% ROBUST (optional) set to 1 or true means to use a slight variation of the
% algorithm that might return duplicates of some intersection points, and
% then remove those duplicates.  The default is true, but since the
% algorithm is slightly slower you can set it to false if you know that
% your curves don't intersect at any segment boundaries.  Also, the robust
% version properly handles parallel and overlapping segments.
%
% The algorithm can return two additional vectors that indicate which
% segment pairs contain intersections and where they are:
%
%   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
% (how far along this segment the intersection is).  For example, if I(k) =
% 45.25 then the intersection lies a quarter of the way between the line
% segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
% the vector J and the segments in (X2,Y2).
%
% You can also get intersections of a curve with itself.  Simply pass in
% only one curve, i.e.,
%
%   [X0,Y0] = intersections(X1,Y1,ROBUST);
%
% where, as before, ROBUST is optional.
%
% Version: 1.10, 25 February 2008
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Modified and Checked in by Shann-Ching Chen, June 12, 2008
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


% Theory of operation:
%
% Given two line segments, L1 and L2,
%
%   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
%   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
%
% we can write four equations with four unknowns and then solve them.  The
% four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
% L1 and L2, t1 is the distance from the starting point of L1 to the
% intersection relative to the length of L1 and t2 is the distance from the
% starting point of L2 to the intersection relative to the length of L2.
%
% So, the four equations are
%
%    (x1(2) - x1(1))*t1 = x0 - x1(1)
%    (x2(2) - x2(1))*t2 = x0 - x2(1)
%    (y1(2) - y1(1))*t1 = y0 - y1(1)
%    (y2(2) - y2(1))*t2 = y0 - y2(1)
%
% Rearranging and writing in matrix form,
%
%  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
%        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
%   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
%        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
%
% Let's call that A*T = B.  We can solve for T with T = A\B.
%
% Once we have our solution we just have to look at t1 and t2 to determine
% whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
% line segments cross and we can include (x0,y0) in the output.
%
% In principle, we have to perform this computation on every pair of line
% segments in the input data.  This can be quite a large number of pairs so
% we will reduce it by doing a simple preliminary check to eliminate line
% segment pairs that could not possibly cross.  The check is to look at the
% smallest enclosing rectangles (with sides parallel to the axes) for each
% line segment pair and see if they overlap.  If they do then we have to
% compute t1 and t2 (via the A\B computation) to see if the line segments
% cross, but if they don't then the line segments cannot cross.  In a
% typical application, this technique will eliminate most of the potential
% line segment pairs.


% Input checks.
narginchk(2,5);

% Adjustments when fewer than five arguments are supplied.
switch nargin
	case 2
		robust = true;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 3
		robust = x2;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 4
		robust = true;
		self_intersect = false;
	case 5
		self_intersect = false;
end

% x1 and y1 must be vectors with same number of points (at least 2).
if length(x1) ~= length(y1)
	error('X1 and Y1 must be equal-length vectors!')
end

if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
        sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1
    x0=[];y0=[];iout=[];jout=[];
    return
end

% x2 and y2 must be vectors with same number of points (at least 2).
if length(x2) ~= length(y2)
	error('X2 and Y2 must be equal-length vectors!!')
end


% Force all inputs to be column vectors.
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
xy1 = [x1 y1];
xy2 = [x2 y2];
dxy1 = diff(xy1);
dxy2 = diff(xy2);

% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2. 
%Converted this from repmat to bsxfun because it is ~20 times faster - HLE
[i,j] = find(bsxfun(@le,min(x1(1:end-1),x1(2:end)),max(x2(1:end-1),x2(2:end)).') & ...
             bsxfun(@ge,max(x1(1:end-1),x1(2:end)),min(x2(1:end-1),x2(2:end)).') & ...
             bsxfun(@le,min(y1(1:end-1),y1(2:end)),max(y2(1:end-1),y2(2:end)).') & ...
             bsxfun(@ge,max(y1(1:end-1),y1(2:end)),min(y2(1:end-1),y2(2:end)).'));

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
% At the same time we can remove redundant combinations of i and j in the
% case of finding intersections of a line with itself.
if self_intersect
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
else
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
end
i(remove) = [];
j(remove) = [];
i = i(:);%Force i and j to be column vectors. Not sure why, but occasionally find is returning them as rows. -HLE
j = j(:);

% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  AA is a 3-D extension of A where we'll use one
% plane at a time.
n = length(i);
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1(i,:).';
AA([2 4],2,:) = dxy2(j,:).';
B = -[x1(i) x2(j) y1(i) y2(j)].';

% Loop through possibilities.  Trap singularity warning and then use
% lastwarn to see if that plane of AA is near singular.  Process any such
% segment pairs to determine if they are colinear (overlap) or merely
% parallel.  That test consists of checking to see if one of the endpoints
% of the curve 2 segment lies on the curve 1 segment.  This is done by
% checking the cross product
%
%   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
%
% If this is close to zero then the segments overlap.

% If the robust option is false then we assume no two segment pairs are
% parallel and just go ahead and do the computation.  If A is ever singular
% a warning will appear.  This is faster and obviously you should use it
% only when you know you will never have overlapping or parallel segment
% pairs.

%Small number for identifying close intersection indices which should be
%integers but are slightly off due to numerical error. 
epsIJ = 1e-6;

if robust
	overlap = false(1,n);
	warning_state = warning('off','MATLAB:singularMatrix');
	% Use try-catch to guarantee original warning state is restored.
	try
		lastwarn('')
		for k = 1:n
			T(:,k) = AA(:,:,k)\B(:,k);
			[unused,last_warn] = lastwarn;
			lastwarn('')
			if strcmp(last_warn,'MATLAB:singularMatrix')
				% Force in_range(k) to be false.
				T(1,k) = NaN;
				% Determine if these segments overlap or are just parallel.
				overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
			end
		end
		warning(warning_state)
    catch err %Replaced lasterr with error catching and rethrowing - HLE
		warning(warning_state)
		rethrow(err)
	end
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	%in_range = T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1;
    %Altered this line to allow for numerical error, which was causing some
    %intersection points to be missed. -HLE
    in_range = T(1,:) >= -epsIJ  & T(2,:) >= -epsIJ & ... %t1 and t2 are ~ greater than 0
               T(1,:) <= (1+epsIJ) & T(2,:) <= (1+epsIJ); %t1 and t2 are ~ less than 1 
    
	% For overlapping segment pairs the algorithm will return an
	% intersection point that is at the center of the overlapping region.
	if any(overlap)
		ia = i(overlap);
		ja = j(overlap);
        T(1,overlap) = .5;%Return the t1 and t2 as .5 to agree with intersections point - HLE                
        T(2,overlap) = .5;
		% set x0 and y0 to middle of overlapping region.
		T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
			min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
		T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
			min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
		selected = in_range | overlap;
	else
		selected = in_range;
	end
	xy0 = T(3:4,selected).';
	
	% Remove duplicate intersection points.
	%[xy0,index] = unique(xy0,'rows'); - Moved this to below the iout &
	%jout calculation - HLE
    index = 1:size(xy0,1);
	x0 = xy0(:,1);
	y0 = xy0(:,2);   
	
	% Compute how far along each line segment the intersections are.
    if nargout > 2
		sel_index = find(selected);
		sel = sel_index(index);
        if ~isempty(sel) %Check if there are any intersections first to avoid error - HLE
            iout = i(sel) + T(1,sel).';
            jout = j(sel) + T(2,sel).';
        else
            iout = [];
            jout = [];
        end

    
        if numel(iout) > 1
            %Remove duplicate intersection points. I Moved this to after
            %computing iout and jout. This allows truly redundant intersections
            %to be separated from instances where there are two distinct
            %intersections at different indices on the polygon which occur at
            %the same point because one polygon either intersects itself or is
            %colinear with itself. Additionally, when two line segments touch
                %at one point, this often is returned as multiple, very close but
            %not identical intersections. This is also handled here so that
            %these are returned as a single intersection. HLE, 7/2011

            %First, simply round any indices which are within epsilon of an
            %integer value. This is simply for convenience of use of the
            %output, as the loop below would still remove any reduntant
            %near-integer indices.
            iout(abs(iout - round(iout)) < epsIJ) = round(iout(abs(iout -round(iout)) < epsIJ));
            jout(abs(jout - round(jout)) < epsIJ) = round(jout(abs(jout -round(jout)) < epsIJ));

            %Loop through the indices, setting those which are closer together
            %than epsilon equal to one another. Repeat this until no indices
            %are this close and non-identical.
            iTooClose = 1;jTooClose =1;
            while ~isempty(iTooClose) || ~isempty(jTooClose)

                [ioutSorted,sortInd] = sort(iout);
                iTooClose = find(diff(ioutSorted) < epsIJ & diff(ioutSorted) > 0);%Find values that are closer than epsilon to each other
                iout(sortInd(iTooClose))  = iout(sortInd(iTooClose+1));%Set these equal to each other, preserving the original order               

                [joutSorted,sortInd] = sort(jout);
                jTooClose = find(diff(joutSorted) < epsIJ & diff(joutSorted) > 0);%Find values that are closer than epsilon to each other
                jout(sortInd(jTooClose))  = jout(sortInd(jTooClose+1));%Set these equal to each other, preserving the original order
            end

            %Now get only the remaining unique intersection points, based on
            %the i and j indices

            %If one of the curves is closed, we need to take into account that
            %the first and last segments are adjacent
            if x1(1) == x1(end) && y1(1) == y1(end)            
                iCheck = iout;
                iCheck(mod(iout,numel(x1))==0) = 1;            
            else
                iCheck = iout;
            end
            if x2(1) == x2(end) && y2(1) == y2(end)
                jCheck = jout;
                jCheck(mod(jout,numel(x2))==0) = 1;
            else
                jCheck = jout;
            end

            [~,index] = unique([iCheck jCheck],'rows');        

            iout = iout(index);
            jout = jout(index);
            x0 = x0(index);
            y0 = y0(index);    

        end
    end
        
else % non-robust option
	for k = 1:n
		[L,U] = lu(AA(:,:,k));
		T(:,k) = U\(L\B(:,k));
	end
	
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1;
	x0 = T(3,in_range).';
	y0 = T(4,in_range).';
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		iout = i(in_range) + T(1,in_range).';
		jout = j(in_range) + T(2,in_range).';
	end
end

if nargout > 2 && ~isempty(iout) && ~isempty(jout) %Again, make sure there were intersections first. - HLE
    [~, sortI] = sort(iout,'ascend');
    iout = iout(sortI); 
    idx = isnan(iout) | isnan(jout);  iout(idx) = [];
    jout = jout(sortI); jout(idx) = [];
    x0 = x0(sortI);     x0(idx) = [];
    y0 = y0(sortI);     y0(idx) = [];
end
% Plot the results (useful for debugging).
%plot(x1,y1,x2,y2,x0,y0,'ok');
