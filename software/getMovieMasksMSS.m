function movieData = getMovieMasksMSS(movieDataOrProcess,varargin)
%THRESHOLDMOVIE applies multi-scale steerable filtering to every frame in input movie
%
% movieData = thresholdMovie(movieData,paramsIn)
%
% Applies manual or automatic thresholding to every frame of the input
% movie and then writes the resulting mask to file as a binary .tif in a
% sub-folder of the movie's analysis directory named "masks"
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
% Output:
%
%   movieData - the updated MovieData object with the thresholding
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The masks are written to the directory specified by the parameter
%   OuptuDirectory, with each channel in a separate sub-directory. They
%   will be stored as binary, bit-packed, .tif files. 
%
%
% Sebastien Besson, Sep 2011
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

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieDataOrProcess', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;


% Get MovieData object and Process
[movieData, segProc] = getOwnerAndProcess(movieDataOrProcess,'MSSSegmentationProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(segProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows'),
    wtBar = waitbar(0,'Initializing...','Name',segProc.getName());
else
end

if ~all(segProc.checkChanNum(p.ChannelIndex))
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

%Read various constants
imDirs  = movieData.getChannelPaths();
nFrames=movieData.nFrames_;

% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    if isempty(p.ProcessIndex)
        inFilePaths{1,i} = imDirs{i};
    else
       inFilePaths{1,i} = movieData.processes_{p.ProcessIndex}.outFilePaths_{1,i}; 
    end
end
segProc.setInFilePaths(inFilePaths);
    
% Set up the output file
outputDir = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outputDir{1,i} = [p.OutputDirectory filesep 'mask_for_channel_' num2str(i)];
    mkClrDir(outputDir{1,i})
end
segProc.setOutFilePaths(outputDir);

%% ---------------Mask calculation ---------------%%% 


disp('Starting applying MSS segmentation...')
%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);

% Anonymous functions for reading input/output
outMask=@(chan,frame) [outputDir{chan} filesep 'mask_' numStr(frame) '.tif'];


logMsg = @(chan) ['Please wait, applying multi-scales filter to channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;
for i=1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
    disp(logMsg(iChan))
    disp(imDirs{iChan});
    disp('Results will be saved under:')
    disp(outputDir{iChan});
    
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iChan)); end
    
    for j=1:nFrames
        % Read image apply mask and save the output
        %Load the current image
        if isempty(p.ProcessIndex)
            currImage = movieData.channels_(iChan).loadImage(j);
        else
            currImage = movieData.processes_{p.ProcessIndex}.loadOutImage(iChan,j);
        end
        
        mask=getCellMaskMSS(double(currImage),'Scales',p.Scales,'FilterOrder',p.FilterOrder);        
        imwrite(mask,outMask(iChan,j));

        % Update the waitbar
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
    
end

% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished segmenting!')
