function smoothedMap = smoothActivityMap(matIn,varargin)
%SMOOTHACTIVITYMAP smooths the input 2D matrix and handles missing observations (NaNs)
% 
% smoothedMap = smoothActivityMap(matIn);
% smoothedMap = smoothActivityMap(matIn,'ParamName1',paramVal1,'ParamName2',paramVal2,...)
% 
%   This function returns a prettier, smoothed version of the input 2D
%   matrix (usually an activity map, e.g. protrusion velocity or
%   fluorescence samples from windows), and handles missing observations
%   (NaNs) by either interpolating to fill these areas or by ignoring them
%   in the smoothing process and leaving these areas blank.
%
%   NOTE: This function is intended to be used only for display/figure
%   making, and the resulting smoothed matrix should NOT be used for
%   analysis.
% 
%   Input:
% 
%       matIn - MxN 2D matrix (activity map) to smooth.
% 
%   Optional Parameters:
% 
%       ('SmoothParam' -> Scalar 0 <= x <= 1). The smoothing parameter to to
%       use to smooth the matrix (this is passed to csaps.m). Must be between
%       0 and 1 where 0 is the smoothest and 1 is least smoothing
%       (interpolating). 
%       Optional. Default is .99
%
%       ('UpSample' -> Positive integer scalar >= 1). Factor by which to
%       upsample. Output matrix will be larger than input matrix by this
%       factor.
%       Optional. Default is 5.
%
%       ('FillNaN' -> True/False) If true, areas with NaN will be filled in
%       with interpolated values. If false, smoothing will ignore these
%       areas. It is generally not a good idea to turn this option on
%       (true), because you are hiding holes in the data, unless these are
%       small holes or outliers (segmentation errors etc) which have been
%       removed. Even if it is off the data will still be smoothed around
%       these holes. 
%       Optional. Default is false.
%
%   Output:
%
%       smoothedMap - u*M x u*N smoothed 2D activity matrix, where u is the
%       upsample parameter.
%
% Hunter Elliott
% 9/2012
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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


ip = inputParser;
ip.addRequired('matIn',@(x)(ismatrix(x) && min(size(x)) > 1));
ip.addParamValue('SmoothParam',.99,@(x)(x >= 0 && x <= 1 && isscalar(x)));
ip.addParamValue('UpSample',5,@(x)(isscalar(x) && x >= 1 && isfinite(x) && isequal(round(x),x)));
ip.addParamValue('FillNaN',false,@(x)(isscalar(x)));
ip.parse(matIn,varargin{:});
p = ip.Results;

[M,N] = size(matIn);
x = 0:N-1;
y = 0:M-1;
xi = linspace(0,N-1,N*p.UpSample);
yi = linspace(0,M-1,M*p.UpSample);

origNanLocs = isnan(matIn);

if nnz(origNanLocs) > 0
    %Interpolate NaNs prior to smoothing, otherwise csaps will remove the
    %entire column for each NaN. 
    [X,Y] = meshgrid(x,y);
    
    %We need to fill in any NaN on the border so all pixels will be within
    %the triangulation used for interpolation. This is obviously sloppy but
    %since this is only for visualization it should be fine.
    medVal = nanmedian(matIn(:));
    matIn(end,isnan(matIn(end,:))) = medVal;
    matIn(1,isnan(matIn(1,:))) = medVal;
    matIn(isnan(matIn(:,1)),1) = medVal;
    matIn(isnan(matIn(:,end)),end) = medVal;
    nanLocs = isnan(matIn);
    
    %Use non-gridded interpolatioin, as the interp2 function doesn't handle
    %NaNs
    intFun = TriScatteredInterp(X(~nanLocs),Y(~nanLocs),matIn(~nanLocs));
    [yNan,xNan] = find(nanLocs);
    intVals = intFun(xNan-1,yNan-1);
    matIn(nanLocs) = intVals;    
    
    %Disable the warning we know the spline function is going to give
    warning('off','SPLINES:CHCKXYWP:NaNs')
    
end

%Do the smoothing on the filled-in matrix. Otherwise the csaps (and other
%spline) function will exclude entire columns which have NaNs in them.
smoothedMap = csaps({y, x},matIn,p.SmoothParam,{yi, xi});

%If we don't want to fill these in, then remove the areas which were NaNs now that
%the matrix has been smoothed.
if ~p.FillNaN && nnz(origNanLocs) > 0   
    %Scale up the original nan matrix so these can be removed from the
    %upsampled activity map. Do it in a slightly conservative way.
    [Xi,Yi] = meshgrid(xi,yi);
    nanLocs = interp2(X,Y,single(origNanLocs),Xi,Yi,'nearest');
    smoothedMap(nanLocs>0) = NaN;    
    
    %Restore warning state
    warning('on','SPLINES:CHCKXYWP:NaNs')
end




