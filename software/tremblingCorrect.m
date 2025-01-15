function tremblingCorrect(movieDataOrProcess, varargin)
% tremblingCorrect wrapper function for TremblingCorrectionProcess.
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

%% ------------------ Input ---------------- %%
ip = inputParser;
ip.addRequired('MD', @(x) isa(x,'MovieData') || isa(x,'Process') && isa(x.getOwner(),'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess, varargin{:});
paramsIn = ip.Results.paramsIn;

%% Registration
% Get MovieData object and Process
[movieData, thisProc] = getOwnerAndProcess(movieDataOrProcess, 'TremblingCorrectionProcess', true);
p = parseProcessParams(thisProc, paramsIn); % If parameters are explicitly given, they should be used
% rather than the one stored in TremblingCorrectionProcess

% Parameters:
currOutputDirectory = p.OutputDirectory;
currRefinementClosureRadius = p.refinementClosureRadius;
currSegProcessIndex = p.SegProcessIndex;

% Sanity Checks
nChan = numel(movieData.channels_);
if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

% precondition / error checking
%Make sure the movie has been segmented
if isempty(currSegProcessIndex)
    iSegProc = movieData.getProcessIndex('MaskProcess',1,true); % nDesired = 1 ; askUser = true
elseif isa(movieData.processes_{currSegProcessIndex},'MaskProcess')
    iSegProc = currSegProcessIndex;
else
    error('The process specified by SegProcessIndex is not a valid MaskProcess! Check input!')
end

if isempty(iSegProc) 
    error('Must create masks before calculating trembling correction! Movie process array has no valid MaskProcess!')
else
   %Check which channels have masks, and use only those that do.
   hasMasks = movieData.processes_{iSegProc}.checkChannelOutput;
   if ~any(hasMasks(p.ChannelIndex))
       error('None of the selected channels have valid masks!');
   end
   if any(~hasMasks(p.ChannelIndex))
        warning('blackWindow:tremblingcorrection:NoMask',...
            'Not all selected channels have masks - using only channels with valid masks!')
        p.ChannelIndex = p.ChannelIndex(hasMasks(p.ChannelIndex));
   end   
end

% logging input paths (bookkeeping)
% input must be some SegProc's output.
inFilePaths = cell(1, numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.processes_{iSegProc}.outFilePaths_{1,i}; 
end
thisProc.setInFilePaths(inFilePaths);

% logging output paths.
dName = 'tremblingCrct_for_channel_';%String for naming the trembling correction directories for each channel
mkClrDir(currOutputDirectory, false);
outFilePaths = cell(1, numel(movieData.channels_));
for iChan = p.ChannelIndex;
    % Create string for current directory
    currDir = [p.OutputDirectory filesep dName num2str(iChan)];
    outFilePaths{1,iChan} = currDir;
    mkClrDir(outFilePaths{1,iChan});
end
thisProc.setOutFilePaths(outFilePaths);


%% Algorithm
% see wrapper2_MSAmultiSeg_imDir_tremblingMaskCorrection.m
for k = p.ChannelIndex
    maskDir = movieData.processes_{iSegProc}.outFilePaths_{k};
    outputDir = outFilePaths{1,k};

    tremblingMaskCorrectionS131(maskDir, outputDir, 'refinementClosureRadius', currRefinementClosureRadius)
end
%%%% end of algorithm

fprintf('\n Finished Trembling Correction! \n')

end