classdef ExternalSegmentationProcess < SegmentationProcess
    % A concrete class for importing segmentation generated by 3d party
    % software
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
    
    methods(Access = public)
        
        function obj = ExternalSegmentationProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            
            super_args{1} = owner;
            super_args{2} = ExternalSegmentationProcess.getName();
            super_args{3} = @importExternalSegmentation;
            if isempty(ip.Results.funParams)
                super_args{4} = ExternalSegmentationProcess.getDefaultParams(...
                    owner, ip.Results.outputDir);
            else
                super_args{4} = ip.Results.funParams;
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end

        function sanityCheck(obj)
            sanityCheck@SegmentationProcess(obj)

            p = obj.getParameters();
            % Test valid channel index matches input data
            validChannels = find(~cellfun(@isempty, p.InputData));
            assert(isequal(validChannels(:), p.ChannelIndex(:)), 'lccb:set:fatal', ...
                'Selected channels do not match input data\n');

            for i = p.ChannelIndex
                if ~exist(p.InputData{i}, 'dir')
                    error('lccb:set:fatal', ...
                        ['The specified mask directory:\n\n ',p.InputData{i}, ...
                        '\n\ndoes not exist. Please double check your channel path.'])
                end
                fileNames = imDir(p.InputData{i},true);
                if isempty(fileNames)
                    error('lccb:set:fatal', ...
                        ['No proper mask files are detected in:\n\n ',p.InputData{i}, ...
                        '\n\nPlease double check your channel path.'])
                end
            end
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'External Segmentation';
        end
        
        function h = GUI()
            h= @externalSegmentationProcessGUI;
        end

        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'externalSegmentation'];
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.InputData = cell(numel(owner.channels_), 1);
        end
    end
end
