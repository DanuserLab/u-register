classdef WindowingPackage < Package
    % The main class of the Windowing package
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
    
    % Sebastien Besson, July 2011
    
    methods
        function obj = WindowingPackage(owner,varargin)
            % Constructor of class WindowingPackage
            
            if nargin == 0
                super_args = {};
            else
                 % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieObject'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'WindowingPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        function parentID = getParent(obj,procID)
            % Use default getParent method
            parentID=getParent@Package(obj,procID);
            
            % Refine dependency relationship between protrusion and
            % windowing processes
            windProcID = find(cellfun(@(x) isa(x,'WindowingProcess'), obj.processes_));
            ProtProcID = find(cellfun(@(x) isa(x,'ProtrusionProcess'), obj.processes_));

            if procID == windProcID
                if isempty(obj.processes_{procID})
                    parentID(parentID==ProtProcID) = [];
                elseif ~strcmp(obj.processes_{procID}.funParams_.MethodName, ...
                        'ProtrusionBased')
                    parentID(parentID==ProtProcID) = [];
                end
            end
        end
    end
    
    methods (Static)
        
        function m = getDependencyMatrix(i,j)
            
            m = [0 0 0 0 0 0 0 0;  %1 GenerateSummationChannelProcess (optional)
                 2 0 0 0 0 0 0 0;  %2 Segmentation
                 0 1 0 0 0 0 0 0;  %3 TremblingCorrectionProcess (optional)
                 0 1 2 0 0 0 0 0;  %4 MaskRefine (optional)
                 0 1 0 2 0 0 0 0;  %5 ProtrusionProcess (optional)
                 0 1 0 2 2 0 0 0;  %6 WindowingProcess
                 0 1 0 2 1 1 0 0;  %7 ProtrusionSamplingProcess
                 0 1 0 2 0 1 0 0;];%8 WindowSamplingProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function name = getName()
            name = 'u-register'; % Updated 2024-05-03. The old name was Windowing.
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = windowingPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            windowingProcConstr = {
                @GenerateSummationChannelProcess,...
                @MultiScaleAutoSegmentationProcess,...
                @TremblingCorrectionProcess,...
                @MaskRefinementProcess,...
                @ProtrusionProcess,...
                @WindowingProcess,...
                @ProtrusionSamplingProcess,...
                @WindowSamplingProcess};
            
            if nargin==0, index=1:numel(windowingProcConstr); end
            procConstr=windowingProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            windowingClasses = {
                'GenerateSummationChannelProcess',...
                'SegmentationProcess',...
                'TremblingCorrectionProcess',...
                'MaskRefinementProcess',...
                'ProtrusionProcess',...
                'WindowingProcess',...
                'ProtrusionSamplingProcess',...
                'WindowSamplingProcess'};
            if nargin==0, index=1:numel(windowingClasses); end
            classes=windowingClasses(index);
        end
        
    end
    
    
end

