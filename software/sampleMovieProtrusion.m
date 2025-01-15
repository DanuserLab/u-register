function movieData = sampleMovieProtrusion(movieData,paramsIn)
%SAMPLEMOVIEPROTRUSION samples the movie protrusion corresponding to each window to produce an activity map
% 
% movieData = sampleMovieProtrusion(movieData)
% movieData = sampleMovieProtrusion(movieData,paramsIn)
% 
% This function samples the protrusion vectors which correspond to each
% window for the input movie. The movie must have been windowed and had
% protrusion vectors calcualted prior to running this function. The result
% is a matrix of protrusion/retraction statistics for each window on the
% edge of the cell and for each frame of the movie. This matrix is often
% referred to as an "activity map."
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
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the samples to. If not
%       input, the samples will be saved to the same directory as the
%       movieData, in a sub-directory called "protrusion_samples"
%
%       ('BatchMode' -> True/False) If true, graphical output and user
%       interaction is supressed (i.e. progress bars, dialog and question
%       boxes etc.)
%
%
% Output:
%
%   movieData - The updated MovieData object, with the parameters and
%   locations of the samples stored in it.
% 
%   Additionally, the protrusion samples for will be written to a
%   sub-directory of the OutputDirectory, as a single .mat file.
% 
% 
% Hunter Elliott
% 8/2010
%
%% --------- Parameters ------------ %%
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

pString = 'protrusion_samples';%Prefix for saving samples to file

showPlots = false; %If figures showing the sampling are displayed. For debugging/testing purposes.

%% ---------- Input ---------------- %%


if nargin < 1 || ~isa(movieData,'MovieData')   
    error('The first input must be a valid MovieData object!');        
end

if nargin < 2 || isempty(paramsIn)
    paramsIn = [];  
end
   
%Check if the movie has been sampled before
iProc = movieData.getProcessIndex('ProtrusionSamplingProcess',1,false);
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(ProtrusionSamplingProcess(movieData,movieData.outputDirectory_));
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%Make sure the movie has been windowed, and find the desired process.
iWinProc = movieData.getProcessIndex('WindowingProcess',1,~p.BatchMode);
if isempty(iWinProc)
    error('The movie could not be sampled, because it has not been windowed yet!')
end

%Make sure that the windows are OK
if ~movieData.processes_{iWinProc}.checkChannelOutput;
    error('No movie channels have valid windows!')
end

%Make sure the protrusion vectors have been calculated
iProtProc = movieData.getProcessIndex('ProtrusionProcess',1,~p.BatchMode);
if isempty(iProtProc) || ~movieData.processes_{iProtProc}.checkChannelOutput
    error('The movie must have protrusion vectors calculated before they can be sampled! Please run getMovieProtrusion.m prior to protrusion sampling!')
end



%% -------- Init ---------- %%

nFrames = movieData.nFrames_;


%Set up and store the output directories for the window samples.
mkClrDir(p.OutputDirectory)

%Load the protrusion vectors
protVecs = movieData.processes_{iProtProc}.loadChannelOutput;
%Separate for readability
protrusion = protVecs.protrusion;
smoothedEdge = protVecs.smoothedEdge;
normals = protVecs.normals;


%Initialize sample array
protSamples.avgNormal = nan(movieData.processes_{iWinProc}.nSliceMax_,nFrames);              
protSamples.stdNormal = protSamples.avgNormal;
protSamples.medNormal = protSamples.avgNormal;
protSamples.minNormal = protSamples.avgNormal;
protSamples.maxNormal = protSamples.avgNormal;
protSamples.avgMagnitude = protSamples.avgNormal;
protSamples.stdMagnitude = protSamples.avgNormal;
protSamples.medMagnitude = protSamples.avgNormal;
protSamples.minMagnitude = protSamples.avgNormal;
protSamples.maxMagnitude = protSamples.avgNormal;
protSamples.avgVector = nan(movieData.processes_{iWinProc}.nSliceMax_,nFrames,2);
protSamples.avgNormalVector = nan(movieData.processes_{iWinProc}.nSliceMax_,nFrames,2);
protSamples.n = protSamples.avgNormal;


if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, sampling protrusion...');
end        


%% -------- Sampling -------- %%


disp('Starting protrusion sampling...')


for iFrame = 1:(nFrames-1)

    %Load the windows for this frame
    windows = movieData.processes_{iWinProc}.loadChannelOutput(iFrame);

    %Insure all normals for this frame have unit length - for some reason
    %Sam's software doesn't return them with unit length...
    magNorm = sqrt(dot(normals{iFrame}',normals{iFrame}')');
    normals{iFrame} = normals{iFrame} ./ repmat(magNorm,1,2);    
    
    %Get the normal component of the protrusion vectors for this frame
    protNorm = dot(protrusion{iFrame}',normals{iFrame}')';
    %And the magnitude of the protrusion vectors
    protMag = sqrt(dot(protrusion{iFrame}',protrusion{iFrame}'))';
    
    for iSlice = 1:numel(windows)
        
        if ~isempty(windows{iSlice}) && ~isempty(windows{iSlice}{1})                        
            
            %Find the protrusion vectors which border this window
            iProts = unique(correspondingIndices(windows{iSlice}{1}{end},smoothedEdge{iFrame}'));
            
            protSamples.n(iSlice,iFrame) = numel(iProts);%Number of prot vecs for this window
            protSamples.avgNormal(iSlice,iFrame) = mean(protNorm(iProts),1);
            protSamples.stdNormal(iSlice,iFrame) = std(protNorm(iProts),[],1);
            protSamples.medNormal(iSlice,iFrame) = median(protNorm(iProts),1);
            protSamples.minNormal(iSlice,iFrame) = min(protNorm(iProts),[],1);
            protSamples.maxNormal(iSlice,iFrame) = max(protNorm(iProts),[],1);
            protSamples.avgMagnitude(iSlice,iFrame) = mean(protMag(iProts),1);
            protSamples.stdMagnitude(iSlice,iFrame) = std(protMag(iProts),[],1);
            protSamples.medMagnitude(iSlice,iFrame) = median(protMag(iProts),1);
            protSamples.minMagnitude(iSlice,iFrame) = min(protMag(iProts),[],1);
            protSamples.maxMagnitude(iSlice,iFrame) = max(protMag(iProts),[],1);
            protSamples.avgVector(iSlice,iFrame,:) = mean(protrusion{iFrame}(iProts,:),1);
            protSamples.avgNormalVector(iSlice,iFrame,:) = mean(normals{iFrame}(iProts,:),1);
            
            if showPlots%For debugging/testing
                
                
                if ~exist('edgePlotted','var') || edgePlotted == false %#ok<UNRCH>
                    
                    plotCols = lines(movieData.processes_{iWinProc}.nSliceMax_);
                    
                    clf;
                    quiver(smoothedEdge{iFrame}(:,1),smoothedEdge{iFrame}(:,2),...
                            protrusion{iFrame}(:,1),protrusion{iFrame}(:,2),0,'k');
                        hold on
                    plot(smoothedEdge{iFrame}(1,1),smoothedEdge{iFrame}(1,2),'rx','MarkerSize',20)
                            
                    edgePlotted = true;

                end
                quiver(smoothedEdge{iFrame}(iProts,1),smoothedEdge{iFrame}(iProts,2),...
                        protrusion{iFrame}(iProts,1),protrusion{iFrame}(iProts,2),0,'Color',plotCols(iSlice,:));

                plotDirection(windows{iSlice}{1},'Color',plotCols(iSlice,:))
                plotDirection(windows{iSlice}{1}(end),'Color',plotCols(iSlice,:))                                    
                
                iMid = min(iProts);
                quiver(smoothedEdge{iFrame}(iMid,1),smoothedEdge{iFrame}(iMid,2),...
                        protSamples.avgVector(iSlice,iFrame,1),...
                    protSamples.avgVector(iSlice,iFrame,2),0,'Color','k','LineWidth',2);
                
                if iSlice == numel(windows)
                    edgePlotted = false;
                end
                
            end                        
        
        end                
        
    end        
    
    if ~p.BatchMode && mod(iFrame,5)
        %Update the waitbar (occasionally to minimize slowdown)
        waitbar(iFrame / nFrames,wtBar)
    end



end


%% ----------- Output ----------- %%


% Create process report
[units, scaling] = movieData.processes_{iProc}.getUnits();
procLog= sprintf([...
    'Protrusion sampling summary\n\n',...
    'Mean protrusion\t: %0.3f %s +/- %0.3f %s\n'...
    'Mean retraction\t: %0.3f %s +/- %0.3f %s\n'],...
    nanmean(protSamples.avgNormal(protSamples.avgNormal>0)) * scaling, units,...
    nanstd(protSamples.avgNormal(protSamples.avgNormal>0)) * scaling, units,...
    nanmean(protSamples.avgNormal(protSamples.avgNormal<0)) * scaling, units,...
    nanstd(protSamples.avgNormal(protSamples.avgNormal<0)) * scaling, units);
disp(procLog);
fid=fopen([p.OutputDirectory filesep 'ProtrusionSamplingSummary.txt'],'wt');
fprintf(fid,procLog);
fclose(fid);

%Write the samples to disk
fPath = [p.OutputDirectory filesep pString '.mat'];
save(fPath,'protSamples');
movieData.processes_{iProc}.setOutFilePath(fPath);

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar);
end

%Update the movie data, save it
movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished protrusion sampling!')



