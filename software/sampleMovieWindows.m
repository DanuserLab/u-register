function movieData = sampleMovieWindows(movieData,paramsIn)
%SAMPLEMOVIEWINDOWS samples the movie images in the sampling windows created with getMovieWindows.m 
% 
% movieData = sampleMovieWindows(movieData)
% movieData = sampleMovieWindows(movieData,paramsIn)
% 
% This function goes through each frame in the movie and samples the images
% in the areas occupied by each sampling window that has been created using
% getMovieWindows.m. Various statistics for the pixels inside each window
% are calculated, and stored in an array the same size as the window array. 
% 
% 
% Input:
% 
%   movieData - A MovieData object describing the movie to be processed, as
%   created by setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       samples to. Samples for different channels will be saved as
%       files in this directory.
%       If not input, the samples will be saved to the same directory as the
%       movieData, in a sub-directory called "window_samples"
%
%       ('ChannelIndex' -> Positive integer scalar or vector) Optional. The
%       integer index of the channels to sample images from. If not input,
%       all available channels will be sampled.
%
%       ('ProcessIndex' -> Positive integer scalar) Optional. This
%       specifies the index of an ImageProcessingProcess or
%       DoubleProcessingProcess in the movieData's process array to use the
%       output images of for sampling. If not specified, the raw images
%       will be sampled.
%
%
%       ('BatchMode' -> True/False)
%       If true, graphical output and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes etc.)
%
%
% Output:
%
%   movieData - The updated MovieData object, with the parameters and
%   locations of the samples stored in it.
% 
%   Additionally, the window samples for each channel will be written to a
%   file in the OutputDirectory, as .mat files, with one file per channel.
% 
% 
% Hunter Elliott
% 7/2010
%
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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

%% ---------- Input ---------------- %%

if nargin < 1 || ~isa(movieData,'MovieData')   
    error('The first input must be a valid MovieData object!');        
end

if nargin < 2 || isempty(paramsIn)
    paramsIn = [];  
end
   
%Check if the movie has been sampled before
iProc = movieData.getProcessIndex('WindowSamplingProcess',1,false);
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(WindowSamplingProcess(movieData,movieData.outputDirectory_));
end
winSampProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(winSampProc,paramsIn);

%Make sure the movie has been windowed, and find the desired process.
iWinProc = movieData.getProcessIndex('WindowingProcess',1,~p.BatchMode);
assert(~isempty(iWinProc),'The movie could not be sampled, because it has not been windowed yet!')
winProc= movieData.processes_{iWinProc};

%Make sure that the windows are okay.
assert(winProc.checkChannelOutput,'The window files for the input movie are not valid!')  

% Check 
stack=dbstack;
if ~any(strcmp('Process.run',{stack(:).name})); winSampProc.run(); return; end

%% -------- Init ---------- %%
if ~iscell(p.ProcessIndex), p.ProcessIndex={p.ProcessIndex}; end
if ~iscell(p.ChannelIndex), p.ChannelIndex={p.ChannelIndex}; end
if ~iscell(p.OutputName), p.OutputName={p.OutputName}; end

nFrames = movieData.nFrames_;
imSize = movieData.imSize_;

nChan = cellfun(@numel,p.ChannelIndex);
nInput = sum(nChan);

%Set up and store the output directories for the window samples.
mkClrDir(p.OutputDirectory)
outFilePaths =cell(numel(p.ProcessIndex),numel(movieData.channels_));
for i=1:numel(p.ProcessIndex)
    iProc = p.ProcessIndex{i};
    if isempty(iProc)
        pString='Raw images - channel ';
    else
        parentOutput = movieData.processes_{iProc}.getDrawableOutput;
        iOutput = strcmp(p.OutputName{i},{parentOutput.var});
        pString=[parentOutput(iOutput).name ' - channel '];
    end
    for j=1:numel(p.ChannelIndex{i})
        iChan = p.ChannelIndex{i}(j);
        outFilePaths{i,iChan} =  [p.OutputDirectory filesep pString num2str(iChan) '.mat'];
    end
end
winSampProc.setOutFilePaths(outFilePaths);

%Initialize sample array
sampledFields = {'avg','std','max','min','med'};
nFields = numel(sampledFields);
allSamples(nInput,1)=struct();
for j = 1:nFields
    [allSamples.(sampledFields{j})] = deal(nan(winProc.nSliceMax_,winProc.nBandMax_,nFrames));
end

if isempty(p.SegProcessIndex)
    %Get the mask information from the windowing process
    iSegProc = movieData.processes_{iWinProc}.funParams_.SegProcessIndex;    
    p.MaskChannelIndex = movieData.processes_{iWinProc}.funParams_.ChannelIndex;
else
    iSegProc = p.SegProcessIndex;    
    if isempty(p.MaskChannelIndex)
        p.MaskChannelIndex = find(movieData.processes_{iSegProc}.checkChannelOutput);        
    end
end

assert(isa(movieData.processes_{iSegProc},'MaskProcess'),'The segmentation process specified by the windowing process is invalid! Please check settings and re-run windowing!')
assert(~isempty(p.MaskChannelIndex),'The specified segmentation process does not have valid masks for the specified channels!')
%Store these in the parameter structure.
p.SegProcessIndex = iSegProc;



%% --------- Sampling --------- %%

if ~p.BatchMode && feature('ShowFigureWindows')
    wtBar = waitbar(0,'Please wait, sampling windows...');
else 
    wtBar = -1;
end  

disp('Starting window sampling...');


for iFrame = 1:nFrames
     
    %Load the windows
    currWin = winProc.loadChannelOutput(iFrame);    
        
    %Load the mask(s) to use first, so we can combine them and use this to
    %mask every channel.
    currMask = true(imSize);
    for j = p.MaskChannelIndex
        currMask = currMask & movieData.processes_{iSegProc}.loadChannelOutput(j,iFrame);
    end    
    
    %Go through each channel and sample it
    stack2sample=zeros([imSize nInput]);

    for i=1:numel(p.ProcessIndex)
        iProc = p.ProcessIndex{i};
        for j=1:numel(p.ChannelIndex{i})
            iChan = p.ChannelIndex{i}(j);
            iInput = sum(nChan(1:i-1))+j;

            if ~isempty(iProc)
                stack2sample(:,:,iInput) = movieData.processes_{iProc}.loadChannelOutput(iChan,iFrame,'output',p.OutputName{i});
            else
                stack2sample(:,:,iInput) = movieData.channels_(iChan).loadImage(iFrame);
            end
        end
    end
    
    currSamples = sampleStackWindows(currWin,stack2sample,currMask);
    assert(numel(currSamples)==nInput);
    
    %Copy these into the whole-movie array
    currSize = size(currSamples(1).(sampledFields{1}));%all field arrays are same size
    for i=1:nInput
        for j = 1:nFields
            allSamples(i).(sampledFields{j})(1:currSize(1),1:currSize(2),iFrame) ...
                = currSamples(i).(sampledFields{j});
        end
    end
       
        
    if ishandle(wtBar) && mod(iFrame,5)
        %Update the waitbar occasionally  to minimize slowdown
        waitbar(iFrame/nFrames,wtBar)
    end
    
end

%% ------- Output ------- %%

if ishandle(wtBar), close(wtBar); end

for i=1:numel(p.ProcessIndex)
    for j=1:numel(p.ChannelIndex{i})
        iChan = p.ChannelIndex{i}(j);
        iInput = sum(nChan(1:i-1))+j;

        samples = allSamples(iInput);  %#ok<NASGU>
        save(outFilePaths{i,iChan},'samples');
    end
end

%Update the movie data, save it
winSampProc.setPara(p);%We've stored additional parameters, so add to the process structure.

disp('Finished sampling!')
