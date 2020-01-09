classdef WindowsDisplay < MovieDataDisplay
    %Concrete class for displaying windows
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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
    properties
        Color='r';
        FaceAlpha=.2;
        showNum=5;
        ButtonDownFcn = [];
    end
    methods
        function obj=WindowsDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
        
        function h=initDraw(obj, data, tag, varargin)
            
            windowArgs = {obj.Color,'FaceAlpha',obj.FaceAlpha};
            h = plotWindows(data, windowArgs, obj.showNum);
            set(h, 'Tag', tag, 'ButtonDownFcn', obj.ButtonDownFcn);
        end
        
        function updateDraw(obj, h, data)
            tag = get(h(1), 'Tag');
            delete(h);
            obj.initDraw(data, tag);
        end
    end
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@ischar;
            params(2).name='FaceAlpha';
            params(2).validator=@isscalar;
            params(3).name='showNum';
            params(3).validator=@isscalar;
            params(4).name='ButtonDownFcn';
            params(4).validator=@(x) isempty(x) || isa(x, 'function_handle');
            
        end
        function f=getDataValidator()
            f=@iscell;
        end
    end
end
