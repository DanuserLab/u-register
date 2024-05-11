function tremblingMaskCorrectionS131(refinedmaskDir, outputDir, varargin)
% tremblingMaskCorrectionS CORRECT trembling of masks that causes 
% bias in correlation analysis mostly in the outermost layer. It simply
% takes moving averages of level functions of masks and then construct new
% masks which can reduce such bias or trembling artifact. Here 'simple'
% means the weight of the moving average is the simplest, [1 3 1]/5.
%
% Usage:
% tremblingMaskCorrectionS131(fullfile(pwd, 'MSASeg_refined_masks'), 'tremblingCorrected_masks131')  
%
% Updates:
%   J Noh, 2021/07/06. Minor edits on footnotes and default values. 
%
% J Noh, 2017/09/22
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
ip.addParameter('show', true);
ip.addParameter('refinementClosureRadius', 3);

parse(ip, varargin{:})
p = ip.Results;


%% Load refined masks

fileReads = dir(refinedmaskDir);
ind = arrayfun(@(x) (x.isdir == 0), fileReads);

fileNames = {fileReads(ind).name};
frmax = numel(fileNames);

currMask = cell(frmax, 1);
imgStack = [];
for fr = 1:frmax
    currMask{fr} = imread(fullfile(refinedmaskDir, fileNames{fr}));
    imgStack = cat(3, imgStack, currMask{fr});
    fprintf(1, '%g ', fr);
end


%% maskrefinement can be done before trembling correction

if (p.refinementClosureRadius > 0) 

    p.MinimumSize = 100;        
    p.ObjectNumber = 1;
    p.FillHoles = 1;

    %%%% you can adjust %%%%%%%%%%%%%%%%%%%%%
    p.ClosureRadius = p.refinementClosureRadius;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    refinedMask = cell(frmax, 1);
    for fr = 1:frmax
        refinedMask{fr} = maskRefinementCoreFunc(currMask{fr}, p);
    end

    %
    currMask = refinedMask;

end


%% level function

DistIm = cell(fr,1);
DistOutside = cell(fr,1);
maskLevel = cell(fr, 1);

for fr = 1:frmax

DistIm{fr} = bwdist(~currMask{fr});

%DistOutside{fr} = (-1) * bwdist(currMask{fr}); %first try -1, 1 at the bdd
DistOutside{fr} = (-1) * bwdist(currMask{fr}) + 1;
DistOutside{fr}(DistOutside{fr} == 1) = 0;

maskLevel{fr} = DistIm{fr} + DistOutside{fr};

end

%% MA, new mask
mask2 = currMask; 

for fr = 2:(frmax-1)
    % weighted average
    % -> changed on 2017/10/04
    %wvec = [1 2 1];  
    wvec = [1 3 1];  
    % wvec = normpdf([-2,0,2], 0, 1);
    w = wvec./sum(wvec);

    tmp = zeros(size(currMask{1}));
    for k = 1:3
        tmp = tmp + w(k) .* maskLevel{fr - 2 + k};
    end
    % thresholding
    mask2{fr} = (tmp >= 0.5);
    mask2{fr} = mat2gray(mask2{fr});    % make output readable from ImageJ
end

    mask2{1} = mat2gray(mask2{1});
    mask2{frmax} = mat2gray(mask2{frmax}); 

%% save
masksOutDir = fullfile(outputDir);
if ~isdir(masksOutDir); mkdir(masksOutDir); end

%%
disp('====')
disp('writing mask files')
for fr = 1:frmax
    %disp('=====')
    %disp(['Frame: ', num2str(fr)])   
    fprintf(1, '%g ', fr);
    %Write the refined mask to file
    imwrite(mask2{fr}, fullfile(masksOutDir, ['tremblingCrctd_', fileNames{fr}]), 'Compression','none'); % fixed issue of ImageJ cannot open compressed mask. - Qiongjing (Jenny) Zou, Jan 2023
end

%%
if p.show
    figure
    for fr=1:frmax
        imshow(mask2{fr})
    end
end

end
