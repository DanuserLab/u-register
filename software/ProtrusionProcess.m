classdef ProtrusionProcess < ImageAnalysisProcess
    %
    % Process Class for calculating protrusion vectors using the
    % getMovieProtrusion.m wrapper function.
    %
    % Hunter Elliott
    % 8/2010
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
        
        function obj = ProtrusionProcess(owner,varargin)
            
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
                super_args{2} = ProtrusionProcess.getName;
                super_args{3} = @getMovieProtrusion;
                
                if isempty(funParams)
                    funParams=ProtrusionProcess.getDefaultParams(owner,outputDir);
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
                    'The file specified as output for the function is invalid!')
            else
                obj.outFilePaths_ = filePath;
            end
            
        end
        
        function status = checkChannelOutput(obj)
            %Overrides the generic function - there is only one set of prot
            %vectors for all channels.
            status = false;
            if exist(obj.outFilePaths_,'file')
                tmp = load(obj.outFilePaths_);
                if isfield(tmp,'protrusion') && isfield(tmp,'normals') ...
                        && isfield(tmp,'smoothedEdge')
                    status = true;
                    
                end
            end
        end
        
        function varargout = loadChannelOutput(obj,varargin)
            
            %Make sure the prot vectors are ok
            if ~checkChannelOutput(obj)
                error('Cannot load the protrusion vectors - they could not be found!')
            end
            
            %             prot = load(obj.outFilePaths_);
            
            outputList = {'protrusion','normals','smoothedEdge'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ProtrusionProcess'));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',{},@(x) all(ismember(x,outputList)));
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
            
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            s = load(obj.outFilePaths_,output{:});
            
            if numel(ip.Results.iFrame)>1,
                if isempty(output)
                    varargout{1}=s;
                else
                    for i=1:numel(output),
                        varargout{i}=s.(output{i});
                    end
                end
            else
                if numel(s.protrusion)>iFrame
                    varargout{1} = [s.smoothedEdge{iFrame} s.smoothedEdge{iFrame}+s.protrusion{iFrame}];
                else
                    varargout{1} = zeros(0,4);
                end
            end
            
        end
        
        
        function h=draw(obj,iFrame,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('iFrame',@isnumeric);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,iFrame,varargin{:})
            
            data=obj.loadChannelOutput(iFrame,'output',ip.Results.output);
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
        
    end
    methods (Static)
        function name =getName()
            name = 'Protrusion';
        end
        function h = GUI()
            h= @protrusionProcessGUI;
        end
        
        function output = getDrawableOutput()
            output(1).name='Protrusion vectors';
            output(1).var={'smoothedEdge','protrusion'};
            output(1).formatData=@(x)[x(:,1:2) x(:,3:4)-x(:,1:2)];
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');
        end
        
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);%Default is to combine masks from all channels
            funParams.SegProcessIndex = [];%No default.
            funParams.OutputDirectory = [outputDir filesep 'protrusion'];
            funParams.DownSample = 50;
            funParams.SplineTolerance = 30;%This is the default in protrusionAnalysis, so I use it also.
            funParams.BatchMode = false;            
        end
    end
end