classdef WindowingProcess < ImageAnalysisProcess
    %WINDOWINGPROCESS is a process for creating sampling windows with getMovieWindows.m
    %
    %
    % Hunter Elliott
    % 7/2010
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
    
    properties (SetAccess = protected, GetAccess = public)
        
        %Window number statistics
        nBandMax_
        nSliceMax_
        
    end
    
    methods (Access = public)
        
        function obj = WindowingProcess(owner,varargin)
            
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
                super_args{2} = WindowingProcess.getName;
                super_args{3} = @getMovieWindows;
                if  isempty(funParams)
                    funParams=WindowingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                
            end
            obj = obj@ImageAnalysisProcess(super_args{:});
      
        end
        
        
        function wins = loadChannelOutput(obj,iFrame,varargin)
            
            if nargin < 2 || isempty(iFrame)
                error('You must specify a frame number to load windows for!');
            elseif round(abs(iFrame)) ~= iFrame || iFrame > obj.owner_.nFrames_
                error('The frame number must be a positive integer less than or equal to the number of frames in the movie!')
            end
            
            winNames = dir([obj.outFilePaths_ filesep '*.mat']);
            
            if numel(winNames) < iFrame
                error('The window file for the specified frame does not exist!')
            end
            
            tmp = load([obj.outFilePaths_ filesep winNames(iFrame).name]);
            fNames = fieldnames(tmp);
            if numel(fNames) ~= 1
                error('Invalid window file !');
            end
            wins = tmp.(fNames{1});
            
        end
        
        function setOutFilePath(obj,filePath)
            %Overloads the method from ImageAnalysisProcess because there
            %is only one set of windows for all channels.
            
            if ~exist(filePath,'dir')
                error('lccb:set:fatal',...
                    'The directory specified for output is invalid!')
            else
                obj.outFilePaths_ = filePath;
                
            end
            
            
        end
        function status = checkChannelOutput(obj)
            %Overloads the generic function - only one set for all channels
            %Make sure the windows exist and are valid
            if ~exist(obj.outFilePaths_,'dir') || ...
                    numel(dir([obj.outFilePaths_ filesep '*.mat'])) ~= obj.owner_.nFrames_;
                status = false;
            else
                status = true;
            end
            
        end
        
        function setWinStats(obj,nSliceMax,nBandMax)
            
            if nargin < 3 || isempty(nBandMax) || ...
                    isempty(nSliceMax)  || ...
                    any([nBandMax nSliceMax]<1)
                error('Invalid window numbers!')
            end
            
            obj.nSliceMax_ = nSliceMax;
            obj.nBandMax_ = nBandMax;
            
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
            name = 'Windowing';
        end
        function h = GUI()
            h= @windowingProcessGUI;
        end
        
        function methods = getMethods()
            methods(1).description = 'Constant number';
            methods(1).name = 'ConstantNumber';
            methods(2).description = 'Constant width';
            methods(2).name = 'ConstantWidth';
            methods(3).description = 'Protrusion based';
            methods(3).name = 'ProtrusionBased';
        end

        function output = getDrawableOutput()
            output(1).name='Windows';
            output(1).var='windows';
            output(1).formatData=[];
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@WindowsDisplay;
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
            funParams.SegProcessIndex = []; %No Default.
            funParams.OutputDirectory = [outputDir  filesep 'windows'];
            funParams.MethodName = 'ConstantNumber';
            funParams.PDEPar = []; %No default.
            funParams.NonLinearSolver = 'off';
            funParams.ParaSize = 10;
            funParams.MeshQuality = []; %Use function default
            funParams.PerpSize = 10;
            funParams.ReInit = Inf;
            funParams.StartPoint = []; %No default
            funParams.StartContour = 2;%Use getMaskWindows default.
            funParams.MinSize = 10; %Minimum number of pixels a mask object must have to be windowed.
            funParams.BatchMode = false;
            funParams.StartPointPropag = true;%Propagates the first window starting point from one frame to the next
        end
    end
end