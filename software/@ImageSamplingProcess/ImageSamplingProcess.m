classdef ImageSamplingProcess < ImageAnalysisProcess
    %Abstract process for sampling image intensities
    %
    % Hunter Elliott
    % 4/2013
    %
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
    
    methods (Access = public)
        
        function obj = ImageSamplingProcess(varargin)                        
            obj = obj@ImageAnalysisProcess(varargin{:});            
        end
        
        
        function samp = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'','avg'};
            nOutput = numel(obj.funParams_.ProcessIndex);
            ip =inputParser;
            ip.addRequired('iChan',@(x) obj.checkChanNum(x));
            ip.addOptional('iOutput',1,@(x) ismember(x,1:nOutput));
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iOutput =ip.Results.iOutput;
            output=ip.Results.output;
            
            assert(obj.checkChannelOutput(iChan,iOutput),'Selected channel and process does not have sampling file!')
            
            s = load(obj.outFilePaths_{iOutput,iChan});
            fNames = fieldnames(s);
            assert(numel(fNames) == 1,'Invalid sample file !');
            samp = s.(fNames{1});
            
            if ~isempty(output), samp=samp.(output); end
        end
        
        function status = checkChannelOutput(obj,varargin)
            
            %Checks if the selected channels have valid output files
            nChanTot = numel(obj.owner_.channels_);
            nOutput = numel(obj.funParams_.ProcessIndex);
            ip = inputParser;
            ip.addOptional('iChan',1:nChanTot,@(x)(all(obj.checkChanNum(x) )));
            ip.addOptional('iOutput',1:nOutput,@(x) ismember(x,1:nOutput));
            ip.parse(varargin{:});
            iChan=ip.Results.iChan;
            iOutput=ip.Results.iOutput;
            
            %Makes sure there's at least one .mat file in the speified
            %directory
            if ~isempty(iOutput)
                status =  cellfun(@(x)logical(exist(x,'file')),obj.outFilePaths_(iOutput,iChan));
            else
                status = cellfun(@(x)logical(exist(x,'file')),obj.outFilePaths_(iChan));
            end
        end
        
        function h=draw(obj,iChan,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('iChan',@isscalar);
            ip.addOptional('iOutput',1,@isscalar);
            ip.KeepUnmatched = true;
            ip.parse(iChan,varargin{:})
            iOutput = ip.Results.iOutput;
            
            data=obj.loadChannelOutput(iChan,iOutput,'output','avg');
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{iOutput,iChan}));
            catch ME
                obj.displayMethod_{iOutput,iChan}=...
                    outputList(iOutput).defaultDisplayMethod(iChan);
            end
            
            % Delegate to the corresponding method
            tag = ['process' num2str(obj.getIndex) '_channel' num2str(iChan) '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput,iChan}.draw(data,tag,drawArgs{:});
        end
        
        function output = getSampledOutput(obj)
            % Construct linear output of time series (to be used by
            % integrator)
            processIndex = obj.funParams_.ProcessIndex;
            channelIndex = obj.funParams_.ChannelIndex;
            outputName = obj.funParams_.OutputName;
            if ~iscell(processIndex), processIndex={processIndex}; end
            if ~iscell(channelIndex), channelIndex={channelIndex}; end
            if ~iscell(outputName), outputName={outputName}; end
            
            nChanPerOutput = cellfun(@numel,channelIndex);
            nOutputTot = sum(nChanPerOutput);
            output(nOutputTot,1)=struct();
            for i=1:numel(processIndex)
                procId = processIndex{i};
                if isempty(procId)
                    processName = '';
                    name='Raw images';
                    
                else
                    processName = class(obj.owner_.processes_{procId});
                    parentOutput = obj.owner_.processes_{procId}.getDrawableOutput;
                    iOutput = strcmp(outputName{i},{parentOutput.var});
                    name=parentOutput(iOutput).name;
                end
                
                for j =1:numel(channelIndex{i})
                    iInput = sum(nChanPerOutput(1:i-1))+j;
                    output(iInput).processName = processName;
                    output(iInput).processIndex = procId;
                    output(iInput).channelIndex = channelIndex{i}(j);
                    output(iInput).var = 'avg';
                    output(iInput).name = [name ' - Channel ' num2str(channelIndex{i}(j))];
                end
            end
        end
        
        function drawableOutput = getDrawableOutput(obj)
            % Build the list of drawable output from 
            processIndex = obj.funParams_.ProcessIndex;
            outputName = obj.funParams_.OutputName;
            if ~iscell(processIndex), processIndex={processIndex}; end
            if ~iscell(outputName), outputName={outputName}; end
            
            % Initialize drawable output array
            nOutput = numel(processIndex);
            drawableOutput(nOutput,1)=struct();
            
            timeInterval = obj.owner_.timeInterval_;
            if ~isempty(timeInterval)
                scaling = [1 timeInterval 1];
                xlabel = 'Time (s)';
            else
                scaling = [1 1 1];
                xlabel = 'Frame number';
            end
            
            for i=1:nOutput
                procId = processIndex{i};
                if isempty(procId)
                    % No procId - sampled raw images 
                    drawableOutput(i).name = 'Raw images';
                else
                    % Read output name from parent process & ouput variable
                    parentOutput = obj.owner_.processes_{procId}.getDrawableOutput;
                    iOutput = strcmp(outputName{i},{parentOutput.var});
                    drawableOutput(i).name=parentOutput(iOutput).name;
                end
                
                %
                drawableOutput(i).var='avg';
                drawableOutput(i).formatData=@(x) permute(x,[1 3 2]);
                drawableOutput(i).type='sampledGraph';
                
                % Get process-defined colormap if applicable
                cmap = @(varargin)jet(2^8);
                if ~isempty(procId) && ismember('getColormap',methods(obj.owner_.processes_{procId}))
                    cmap=@(x,i)obj.owner_.processes_{procId}.getColormap(iOutput,obj.getIntensityLimits(x,i));
                end
                
                % Get process-defined units if applicable
                units =@(varargin) '';
                if ~isempty(procId) && ismember('getUnits',methods(obj.owner_.processes_{procId}))
                    units=@(i)obj.owner_.processes_{procId}.getUnits(i);
                end
                
                drawableOutput(i).defaultDisplayMethod = @(x) ScalarMapDisplay(...
                    'Colormap', cmap(x,i), 'Units', units(i),...
                    'CLim', obj.getIntensityLimits(x,i), 'Scaling', scaling,...
                    'Labels', {xlabel, 'Window depth', 'Window number'});
            end
        end
        
        

        
        function limits = getIntensityLimits(obj,iChan,iOutput)
            data=obj.loadChannelOutput(iChan,iOutput,'output','avg');
            limits=[min(data(:)) max(data(:))];
        end
    end

end