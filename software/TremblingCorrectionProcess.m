classdef TremblingCorrectionProcess < MaskProcessingProcess
    % A concrete process Trembling Correction
    % CORRECT trembling of masks that causes bias in correlation analysis mostly in the outermost layer.
    % see tremblingMaskCorrectionS131.m, maskRefinementCoreFunc.m, tremblingCorrect.m
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
    
    % Qiongjing (Jenny) Zou, Aug 2021
    
    % Bug fixed, change super class to MaskProcessingProcess, so user will
    % be asked in MaskRefinementProcess whether to choose this proc as seg
    % process. -- Qiongjing (Jenny) Zou, Dec 2022
    
    methods
        function obj = TremblingCorrectionProcess(owner,varargin)
            
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
                super_args{2} = TremblingCorrectionProcess.getName;
                super_args{3} = @tremblingCorrect;
                if isempty(funParams)
                    funParams=TremblingCorrectionProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@MaskProcessingProcess(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'Trembling Correction';
        end
        function h = GUI()
            h= @TremblingCorrectionProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x,'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir = ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'TremblingCorrection'];
            funParams.refinementClosureRadius = 5;
            funParams.SegProcessIndex = [];%No default.
        end
    end
end