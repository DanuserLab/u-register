function intPer = apportionIntegers(fracs,total)
%APPORTIONINTEGERS apportions integer values with given fraction and total
% 
% intPer = apportionIntegers(frac,total)
% 
% This function apportions (sub-divides) the input total into integer
% sub-totals so that each sub-total has approximately the specified
% fraction of the total, and the sum of these sub-totals equals the input
% total. It uses Hamilton's method to do so. This problem is impossible to
% do without producing some potentially undesireable results (see Balinkski
% & Young, 1982), but this method will always maintain the input total.
% 
% Input:
%   
%   frac - An Mx1 or 1xM vector containing the fraction of the total to
%   apportion to each of the M output integers. This vector should sum to
%   1. Note that as stated above, these fractions will rarely be perfectly
%   satisfied in order to maintain the total.
% 
%   total - An integer specifying the desired sum of the M output integers.
% 
% 
% Output:
% 
%   intPer - The integer values apportioned based on each of the M input
%   fractions, having a sum equal to the input total
%
% Hunter Elliott
% 6/2011
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


if nargin < 2 || isempty(fracs) || isempty(total)
    error('Must input both a list of fractions and a total!')
end

if any(fracs < 0) || numel(fracs) < 2 || ...
        abs(sum(fracs) - 1) > 10*eps(class(fracs)) %Give some leeway in the sum for numerical error
    error('The fracs input must be a positive vector of length greater than 1 which sums to 1');
end

if round(abs(total)) ~= total || numel(total) ~=1
    error('The input total must be a positive integer scalar!')
end

%Get the non-integer fraction each integer should get ("standard quota")
standQuota = fracs .* total;

%Initialy assign the lower quota to each.
lowerQuota = floor(standQuota);
intPer = lowerQuota;

%Determine how many still need to be apportioned
nRem = total - sum(lowerQuota);

%Assign the remaining integers to each output based on the remainder which
%was rounded away when assigning the lower quota
[~,sortRem] = sort(standQuota - lowerQuota,'descend');
intPer(sortRem(1:nRem)) = intPer(sortRem(1:nRem)) + 1;
    




