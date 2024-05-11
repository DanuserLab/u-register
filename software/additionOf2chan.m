function additionOf2chan(movieDataOrProcess, varargin)
% additionOf2chan wrapper function for GenerateSummationChannelProcess.
%
% INPUT
% movieDataOrProcess - either a MovieData (legacy)
%                      or a Process (new as of July 2016)
%
% param - (optional) A struct describing the parameters, overrides the
%                    parameters stored in the process (as of Aug 2016)
%
% OUTPUT
% none (saved to p.OutputDirectory)
%
% Changes
% As of July 2016, the first argument could also be a Process. Use
% getOwnerAndProcess to simplify compatability.
%
% As of August 2016, the standard second argument should be the parameter
% structure
%
% Qiongjing (Jenny) Zou, Aug 2021
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

%% ------------------ Input ---------------- %%
ip = inputParser;
ip.addRequired('MD', @(x) isa(x,'MovieData') || isa(x,'Process') && isa(x.getOwner(),'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess, varargin{:});
paramsIn = ip.Results.paramsIn;

%% Registration
% Get MovieData object and Process
[movieData, thisProc] = getOwnerAndProcess(movieDataOrProcess, 'GenerateSummationChannelProcess', true);
p = parseProcessParams(thisProc, paramsIn); % If parameters are explicitly given, they should be used
% rather than the one stored in GenerateSummationChannelProcess

% Parameters:
currOutputDirectory = p.OutputDirectory;

% Sanity Checks
nChan = numel(movieData.channels_);
if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

% precondition / error checking
if numel(movieData.channels_) < 2
    error('Need at least two channels for this Generate Summation Channel process!')
end

% logging input paths (bookkeeping)
inFilePaths = cell(1, numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
thisProc.setInFilePaths(inFilePaths);

% logging output paths.
mkClrDir(currOutputDirectory, false);
outFilePaths = cell(1, numel(movieData.channels_));
for iChan = p.ChannelIndex;
    % output dir for both chan are the same
    outFilePaths{1,iChan} = currOutputDirectory;
    mkClrDir(outFilePaths{1,iChan});
end
thisProc.setOutFilePaths(outFilePaths);


%% Algorithm
% see additionOf2chanImTiffs.m
% Edit to make it work for all MD.Reader, such as BioFormatsReader. Before, 
% the algorithm only works for TiffSeriesReader. - March 2022

nCh = numel(p.ChannelIndex);
% imgFolderPaths = inFilePaths;
outputDir = currOutputDirectory;
addedChName = 'SummationImg';

I = cell(nCh, 1);
for k = 1:nCh
%     fileReads = dir(imgFolderPaths{1,k});
%     ind = arrayfun(@(x) (x.isdir == 0), fileReads);
    
%     fileNames = {fileReads(ind).name};
    fileNamesF = movieData.getImageFileNames(k); 
    fileNames = fileNamesF{1}; % this is the way to get fileNames which works for diff kinds images, including tiffs, bioFormats images
    frmax = numel(fileNames);
    
%     FileTif= fullfile(imgFolderPaths{1,k}, fileNames{1});
%     InfoImage=imfinfo(FileTif);
%     mImage=InfoImage(1).Width
%     nImage=InfoImage(1).Height
    mImage = movieData.imSize_(2)
    nImage = movieData.imSize_(1)
    
    %I{k} = zeros(nImage, mImage, frmax);
    tmpIm = zeros(nImage, mImage, frmax);
    parfor fr = 1:frmax
%         tmpIm(:,:,fr) = imread(fullfile(imgFolderPaths{1,k}, fileNames{fr}));
        tmpIm(:,:,fr) = movieData.channels_(k).loadImage(fr); % this is the way to read image for all MD.Reader; this works when input is always raw images.
        %imgStack = cat(3, imgStack, I{fr});
        %fprintf(1, '%g ', fr);
        disp(fr)
    end
    I{k} = tmpIm;
end

% Addition of 2channels
I{3} = uint16(zeros(nImage, mImage, frmax));
for fr = 1:frmax
    I{3}(:,:,fr) = I{1}(:,:,fr) + I{2}(:,:,fr);
end

if ~isdir(outputDir); mkdir(outputDir); end
for fr = 1:frmax
    outfname = fullfile(outputDir, [addedChName, '_', sprintf('%04d', fr), '.tif']);
    imwrite(I{3}(:,:,fr), outfname);
    fprintf(1, '%g ', fr);
end

%%%% end of algorithm

fprintf('\n Finished Generate Summation Channel! \n')

end