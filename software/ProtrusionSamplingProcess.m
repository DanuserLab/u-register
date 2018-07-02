classdef ProtrusionSamplingProcess < ImageAnalysisProcess
    %
    % Process Class for sampling the protrusion vectors which correspond with
    % each window for a given movie.
    %
    % Hunter Elliott
    % 1/2011
    %
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
    
    methods (Access = public)
        
        function obj = ProtrusionSamplingProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = ProtrusionSamplingProcess.getName;
                super_args{3} = @sampleMovieProtrusion;                
                if isempty(funParams)
                    funParams=ProtrusionSamplingProcess.getDefaultParams(owner,outputDir);
                end                
                super_args{4} = funParams;
                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        function setOutFilePath(obj,filePath)
            %Overloads the method from ImageAnalysisProcess because there
            %is only one set of vectors for all channels, which is stored
            %as a single file
            
            if ~exist(filePath,'file')
                error('lccb:set:fatal',...
                    'The file specified specified as output for the function is invalid!')
            else
                obj.outFilePaths_ = {filePath};
            end
            
        end
        
        function status = checkChannelOutput(obj)
            %Overrides the generic function - there is only one set of prot
            %samples for all channels.
            status = logical(exist(obj.outFilePaths_{1},'file'));
        end
        
        function prot = loadChannelOutput(obj,varargin)
            
            %Make sure the prot samples are ok
            if ~checkChannelOutput(obj)
                error('Cannot load the protrusion samples - they could not be found!')
            end
            
            outputList = {'','avgNormal'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ProtrusionSamplingProcess'));
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
            ip.parse(obj,varargin{:})
            output = ip.Results.output;
            
            prot = load(obj.outFilePaths_{1});
            fn = fieldnames(prot);
            if numel(fn) > 1, error('Invalid protrusion sample file!'); end
            prot = prot.(fn{1});
            if ~isempty(output), prot=prot.(output); end
            
        end
        
        function h=draw(obj,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,varargin{:})
            
            data=obj.loadChannelOutput('output',ip.Results.output);
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{iOutput}));
            catch ME
                obj.displayMethod_{iOutput}=...
                    outputList(iOutput).defaultDisplayMethod();
            end
            
            % Delegate to the corresponding method
            tag = ['process' num2str(obj.getIndex) '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
        end
        
        function output = getSampledOutput(obj)
            output.processName = 'Protrusion';
            output.channelIndex = obj.owner_.getProcessIndex('ProtrusionProcess',1,~obj.funParams_.BatchMode);
            output.var = 'avgNormal';
            output.name = 'Protrusion';
        end
        function output = getDrawableOutput(obj)
            [units, scaling timeInterval] = obj.getUnits();
            output(1).name='Protrusion map';
            output(1).var='avgNormal';
            output(1).formatData=[];
            output(1).type='movieGraph';
            unitsLabel = ['Edge velocity (' units ') '];
            if ~isempty(timeInterval)
                yxScaling = [1 timeInterval scaling];
                xlabel = 'Time (s)';
            else
                yxScaling = [1 1 1];
                xlabel = 'Frame number';
            end
            output(1).defaultDisplayMethod = @(x)ScalarMapDisplay('Colormap',jet(2^8),...
                'Units', unitsLabel, 'Labels', {xlabel, 'Window number'},...
                'Scaling', yxScaling);
        end

        function [units scaling timeInterval] = getUnits(obj)
            pixelSize = obj.owner_.pixelSize_;
            timeInterval = obj.owner_.timeInterval_;
            if isempty(pixelSize) || isempty(timeInterval),
                units = 'pixels/frame';
                scaling = 1;
            else
                units = 'nm/s';
                scaling = pixelSize/timeInterval;                
            end
        end
    end
    methods (Static)
        function name =getName()
            name = 'Protrusion Sampling';
        end
        function name =GUI()
            name =@noSettingsProcessGUI;
        end

        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'protrusion_samples'];
            funParams.BatchMode = false;
        end

    end
end