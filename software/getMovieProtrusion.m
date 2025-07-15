function movieData = getMovieProtrusion(movieData,paramsIn)
%GETMOVIEPROTRUSION calculates protrusion vectors for the input movie
%
% movieData = getMovieProtrusion(movieData)
%
% movieData = getMovieProtrusion(movieData,paramIn)
% 
% This function is a wrapper for Sam's protrusion calculation function
% prSamProtrusion.m. It allows the protrusion vectors to be calculated
% based on a MovieData object. This enables batch processing and command
% line use.
% 
% In short, the function calculates vectors connecting the edge of the mask
% in frame n with the edge of the mask in frame n+1. Each mask must contain
% only ONE object.
% 
%
% Input:
% 
%   movieData - A MovieData object describing the movie to calculate
%   protrusion vectors for.
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the protrusion vectors to.
%       If not input, the vectors will be saved to the same directory as
%       the movieData, in a sub-directory called "protrusion"
%
%       ('ChannelIndex' -> Positive integer scalar or vector) Optional. The
%       integer index of the channel(s) to use masks from. If not input,
%       all channels with masks will be used. If multiple channels are
%       selected, the intersection of the masks from each channel is used.
% 
%       ('SegProcessIndex' -> Positive integer scalar) This specifies the
%       index of the segmentation process to use masks from forn
%       calculating protrusion. This is only necessary if more than one
%       segmentatin process exists. If not specified, and more than one
%       segmentation process exists, the user will be asked, unless batch
%       mode is enabled, in which case an error will be generated.
%
%       ('DownSample' -> Positive integer scalar) Optional. The
%       downsampling parameter to use for protrusion calculation. Edge
%       segments with more points than this will be downsampled to this
%       number of points. Lower numbers will speed calculation but decrease
%       quality of resulting vectors. Default is 50.
%
%       ('SplineTolerance' -> Positive scalar) Optional. The tolerance
%       value to use when fitting a spline to mask object outlines. Default
%       from protrusionAnalysis is 30. Larger numbers will smooth the mask
%       edge more, while 0 will use an interpolating spline.
%
%       ('BatchMode' -> True/False)
%       If true, graphical output and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes etc.)
%
%
% Output:
%
%   movieData - the updated MovieData object with the
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The protrusion vectors are written to the directory specified by the
%   parameter OuptuDirectory, they are stored as a single .mat file.
%
% 
% Hunter Elliott
% Re-written 8/2010
%
% Daehwan Kim (infphilo@gmail.com)
% 12/2017
% Parallelized protrusion detection
%
% Made it work for cropped movie previously done in Biosensors Package.
% Qiongjing (Jenny) Zou, Jan 2023 & April 2023
%% ------- Parameters ------ %%
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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

fName = 'protrusion_vectors'; %String for naming file with results
samBatch = true; %Wether to run Sam's function in batch mode or not.

%% ----- Input ------ %%

if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end

%Get the indices of any previous mask refinement processes from this function                                                                              
iProc = movieData.getProcessIndex('ProtrusionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(ProtrusionProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);


%Make sure the move has been segmented
if isempty(p.SegProcessIndex)
    iSegProc = movieData.getProcessIndex('MaskProcess',1,~p.BatchMode);
elseif isa(movieData.processes_{p.SegProcessIndex},'MaskProcess')
    iSegProc = p.SegProcessIndex;
else
    error('The process specified by SegProcessIndex is not a valid MaskProcess! Check input!')
end

if isempty(iSegProc) 
    error('Must create masks before calculating protrusion! Movie process array has no valid MaskProcess!')
else
   %Check which channels have masks, and use only those that do.
   hasMasks = movieData.processes_{iSegProc}.checkChannelOutput;
   if ~any(hasMasks(p.ChannelIndex))
       error('None of the selected channels have valid masks!');
   end
   if any(~hasMasks(p.ChannelIndex))
        warning('blackWindow:protrusion:NoMask',...
            'Not all selected channels have masks - using only channels with valid masks!')
        p.ChannelIndex = p.ChannelIndex(hasMasks(p.ChannelIndex));
   end   
end

%% ---- Init ----- %%

%Get mask directories and names
maskDirs = movieData.processes_{iSegProc}.outFilePaths_(p.ChannelIndex);
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.ChannelIndex);

%Set up output directory
mkClrDir(p.OutputDirectory);

nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
% imSize = movieData.imSize_; % does not work for cropped movie

%Set up the inputs for sam's function that are common to all frames
samIn.batch_processing = samBatch;
samIn.TOLERANCE = p.SplineTolerance;
samIn.dl_rate = p.DownSample;

%Disable nearly-singular warning which is produced by intersections.m when
%used by prSamProtrusion.m This warning occurs when attempting to find the
%point of intersection between two line segments. In our case, we only need
%to know which line segements intersect and don't care about the accuracy
%of the intersection point, so it is okay to ignore this warning.
warning('off','MATLAB:nearlySingularMatrix');

%Disable the annoying optimization algorithm options that result from a
%faulty version check in sam's software
if datenum(version('-date')) < datenum('August 13, 2011')
    warning('off','optim:fmincon:NLPAlgLargeScaleConflict');
else
    %For matlab2011b (and presumably future versions) the warning identifier
    %has changed, so we need to use a different command to 
    warning('off','optimlib:fmincon:NLPAlgLargeScaleConflict')
end


%% ----- Protrusion Calculation ---- %%

if ~p.BatchMode
    wtBar = parfor_progressbar(nFrames - 1, 'Please wait, calculating protrusion vectors ...');
end        
disp('Starting protrusion calculation...');

% prevMask = true(imSize); % does not work for cropped movie
isPrevClosed = true; 
protrusion = cell(nFrames-1,1);
normals = cell(nFrames-1,1);
smoothedEdge = cell(nFrames,1);

% Make this process work for previous cropped movie, which is done in CropShadeCorrectROIProcess in BiosensorsPackage. - Jan 2023
inDir = movieData.processes_{p.SegProcessIndex}.outFilePaths_{1,p.ChannelIndex(1)};
dinfo = dir(inDir);
imInfo = imfinfo(fullfile(dinfo(end).folder,dinfo(end).name));
n = imInfo.Height;
m = imInfo.Width;
imSize = [n m];
prevMask = true(n, m);
if isequal(n, movieData.imSize_(1)) && isequal(m, movieData.imSize_(2))
    roiMask=movieData.getROIMask;
else
    roiMask = true([n, m, movieData.nFrames_]);
end

parfor iFrame = 2:nFrames
    [prevMask, isPrevClosed, prevOutline, isPrevClockWise] = ...
        getOutlineInfo(roiMask(:,:,iFrame - 1), iFrame - 1, nChan, maskDirs, maskNames, imSize);
    [currMask, isCurrClosed, currOutline, isClockWise] = ...
        getOutlineInfo(roiMask(:,:,iFrame), iFrame, nChan, maskDirs, maskNames, imSize);
    
    if isCurrClosed ~= isPrevClosed
        error('The function prSamProtrusion.m requries that if the mask object is touching the image boundary in one frame, it must also touch the boundary in every frame! Check masks and re-run calculation!')
    end
    
    %Set up the inputs for sam's function
    samIn = struct();
    samIn.batch_processing = samBatch;
    samIn.TOLERANCE = p.SplineTolerance;
    samIn.dl_rate = p.DownSample;
    
    %Set up input for sam's function
    samIn.pixel_list_last = prevOutline;
    samIn.pixel_list = currOutline;
    samIn.ISCLOSE = isCurrClosed;
    samIn.maskTM1 = prevMask;
    samIn.time_idx = iFrame;
    
    %----call sam's protrusion calc function----%
    samOut = prSamProtrusion(samIn);
        
    %Parse output.
    protrusion{iFrame-1} = samOut.translate_output;
    normals{iFrame-1} = samOut.xy_normal;
    smoothedEdge{iFrame-1} = samOut.pixel_tm1_output;
    if iFrame == nFrames
        smoothedEdge_t{iFrame} = samOut.pixel_t_output;
    end
            
    if ~p.BatchMode        
        wtBar.iterate(1);
    end
end
smoothedEdge{nFrames} = smoothedEdge_t{nFrames};

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


%% ------ Output ----- %%

%Save the protrusion vectors and smoothed cell edge to a file
save([p.OutputDirectory filesep fName '.mat'],...
            'protrusion','normals','smoothedEdge');

movieData.processes_{iProc}.setOutFilePath(...
                                [p.OutputDirectory filesep fName '.mat']);
movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished protrusion calculation!')
end


function [mask, isClosed, outline, isClockWise] = getOutlineInfo(mask, iFrame, nChan, maskDirs, maskNames, imSize)
%Load and combine masks from all channels
for iChan = 1:nChan
    mask = mask & imread([maskDirs{iChan} filesep maskNames{iChan}{iFrame}]);
end

isClosed = ~any([mask(1,:), mask(end,:) mask(:,1)' mask(:,end)']);

%Separate mask objects - prSamProtrusion only handles one object
CC = bwconncomp(mask);

%Make sure there isn't more than one object
if CC.NumObjects > 1
    if nChan == 1
        %Because small objects can be created by intersecting masks, we
        %only warn the user if only one mask was used.
        warning('blackwindow:protrusion:TooManyObjects',...
            ['The mask for frame ' num2str(iFrame) ...
            ' contains more than 1 object - protrusion vectors are only calculated for the largest object.'])
    end
    %We just take the largest of these objects for prot vec calc
    [~,iBiggest] = max(cellfun(@(x)(numel(x)),CC.PixelIdxList));
    mask = false(imSize); % imSize was not a input variable, added 2023-4-13
    mask(CC.PixelIdxList{iBiggest}) = true;
    
end

%Get the outline of the object in this mask. We use contourc instead of
%bwboundaries for 2 reasons: It returns its results in matrix
%coordinates, and the resulting outline encloses the border pixels
%instead of running through their centers. This better agrees with the
%windows, as the windows are designed to enclose the entire mask.
outline = contourc(double(mask),[0 0]);
outline = separateContours(outline);%Post-processing of contourc output
outline = cleanUpContours(outline);
outline = outline{1}';%We know we only have one object...

%Make sure the outline is correctly oriented
if ~isClosed
    %Close the curve before checking handedness
    closedOutline = closeContours({outline'},bwdist(~mask));
    isClockWise = isCurveClockwise(closedOutline{1});
else
    isClockWise = isCurveClockwise(outline);
end

if ~isClockWise
    %Sam requires the curves run in the same direction
    outline = outline(end:-1:1,:);
end
end

