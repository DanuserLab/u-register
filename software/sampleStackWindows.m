function samples = sampleStackWindows(windows,stack,mask)
%SAMPLEIMAGEWINDOWS samples the stack intensity in each of the input windows
%
% samples = sampleStackWindows(windows,image,mask);
%
% This function calculates and returns various statistics regarding the
% pixel values in the input image within each of the input window polygons.
% The input windows should be of the form created by getMaskWindows.m or
% getMovieWindows.m
% 
% Input:
% 
%   windows - Cell-array containing window polygons, as produced by
%   getMaskWindows.m or getMovieWindows.m
% 
%   stack - 3D matrix. The stack to sample. Each slice must be the same size as the
%   mask which was used to create the window array.
% 
%   mask - A 2D logical matrix of the same size as that used to create the
%   windows. Pixel values outside of this mask will treated as NaNs and
%   ignored in calculating image statistics. If you don't want to mask
%   additionally, just use the mask you used to create the windows.
% 
% Output:
% 
%   samples - Structure with fields containing the statistics of the image
%   pixels within each window. To ease analysis, instead of a cell array
%   like the windows, each field is a 2D MxN matrix, where N is the number
%   of window columns ("slices") and M is the maximum number of window rows
%   ("bands") in any column. Therefore, any elements of this matrix which
%   do not have a corresponding window will have NaN in their statistics.
%   Also, any windows which do not cover any masked pixels in the image
%   will have NaN in their statistics.
%
%     Details:
%
%       samples.avg - MxN matrix of average of pixel values in each window.
%       samples.med - MxN matrix of median of pixel values in each window.
%       samples.std - MxN matrix of standard deviation of pixel values in each window.
%       samples.min - MxN matrix of minimum of pixel values in each window.
%       samples.max - MxN matrix of maximum of pixel values in each window.
%       samples.n   - MxN matrix of number of pixels in each window.
%       
%
% Hunter Elliott
% 3/2011
%
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


% NOTE: The note below is old :). Now we allow different sets of mask - one
% for defining the cell outline, and another which may be based on e.g SNR
% in the FRET signal. Using two pairs allows us to track morphodynamics
% with an accurate mask (potentially from another image channel) while
% still excluding cell regions with SNR lowe enough to have unreliable FRET
% ratios / fluorescence signals. When this isn't the case, see the note
% below!
% NOTE: We need to use the masks again here during sampling for several
% reasons. One, neither inpolygon.m nor poly2mask.m produce results which
% agree completely with the output of contourc.m or bwboundaries.m. That
% is, creating a polygon that surrounds a mask using either contourc.m or
% bwboundaries.m and then using poly2mask or inpolygon will NOT re-create
% the same mask. Second, the protrusion calculation functions use a
% smoothed version of the cell edge, which does not faithfully reproduce
% the outline of the mask. Therefore, any protrusion-based window
% propagation methods will no longer give windows which perfectly cover the
% mask. Third, if other window propagation methods are used, the windows
% may not completely agree with the masks in frames after the
% initialization/re-initialization frames. Therefore, as stupid as it
% sounds, we need to re-mask the images before we sample the window
% intensities, to avoid sampling un-masked values.
%   -Hunter


% Input check
ip=inputParser;
ip.addRequired('windows');
ip.addRequired('stack');
ip.addRequired('mask',@(x) islogical(x) && ndims(mask) == 2);
ip.parse(windows,stack,mask);

assert(size(stack,1)==size(mask,1) && size(stack,2) == size(mask,2),...
    'The image and mask must have the same size along the first 2 dimensions!');

%% ------------------- Init ------------------------- %%

%Convert to double to avoid rounding-error if integer-valued images are
%input.
stack = double(stack);

imSize = size(stack);

nSlices = numel(windows);
nBands = cellfun(@numel,windows);

%Initialize sample array
sampledFields = {'avg','med','std','min','max','n'};
nFields = numel(sampledFields);
samples(size(stack,3),1)=struct();
for j = 1:nFields
    [samples.(sampledFields{j})] = deal(nan(nSlices,max(nBands)));
end


%% -------------------- Sampling -------------------- %%

for iSlice = 1:nSlices    
    for iBand = 1:nBands(iSlice)
        
       %Make sure it's not an empty placeholder window 
       if ~isempty(windows{iSlice}{iBand})
           
            %Create mask for this window. Windows are in matrix coord.
            windowsPoly = [windows{iSlice}{iBand}{:}];

            winMask = poly2mask(windowsPoly(1,:),...
                                windowsPoly(2,:),...
                                imSize(1),imSize(2)); %Faster to create mask that only covers this window?? Faster to use inpolygon?                                    

            %Remove any un-masked pixels in any of the mask
            %channel(s). See note above as to why this is
            %necessary.
            winMask = winMask & mask;
            mask(winMask) = false; %Remove these pixels from the mask. If there's something wrong with the propagation method, we don't want ot over-sample pixels.

            if nnz(winMask)>0            
                for nz=1:size(stack,3)
                    im=stack(:,:,nz);
                    %Sample the image with this mask
                    currSamp = im(winMask(:));
                    
                    %Calculate statistics of this sample
                    samples(nz).avg(iSlice,iBand) = nanmean(currSamp);
                    samples(nz).med(iSlice,iBand) = nanmedian(currSamp);
                    samples(nz).std(iSlice,iBand) = nanstd(currSamp);
                    samples(nz).min(iSlice,iBand) = nanmin(currSamp);
                    samples(nz).max(iSlice,iBand) = nanmax(currSamp);
                    samples(nz).n(iSlice,iBand) = numel(currSamp);
                end
            end            
       end
        
    end    
end




