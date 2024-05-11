function varargout = GenerateSummationChannelProcessGUI(varargin)
%GENERATESUMMATIONCHANNELPROCESSGUI MATLAB code file for GenerateSummationChannelProcessGUI.fig
%      GENERATESUMMATIONCHANNELPROCESSGUI, by itself, creates a new GENERATESUMMATIONCHANNELPROCESSGUI or raises the existing
%      singleton*.
%
%      H = GENERATESUMMATIONCHANNELPROCESSGUI returns the handle to a new GENERATESUMMATIONCHANNELPROCESSGUI or the handle to
%      the existing singleton*.
%
%      GENERATESUMMATIONCHANNELPROCESSGUI('Property','Value',...) creates a new GENERATESUMMATIONCHANNELPROCESSGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to GenerateSummationChannelProcessGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GENERATESUMMATIONCHANNELPROCESSGUI('CALLBACK') and GENERATESUMMATIONCHANNELPROCESSGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GENERATESUMMATIONCHANNELPROCESSGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

% Edit the above text to modify the response to help GenerateSummationChannelProcessGUI

% Last Modified by GUIDE v2.5 11-Aug-2021 13:46:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GenerateSummationChannelProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GenerateSummationChannelProcessGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GenerateSummationChannelProcessGUI is made visible.
function GenerateSummationChannelProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Parameter setup
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
funParams = userData.crtProc.funParams_;


% set GUI with Parameters
% no parameters for this GUI.

% Update user data and GUI data
handles.output = hObject;
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = GenerateSummationChannelProcessGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)


%  Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String')) || ~isequal(numel(get(handles.listbox_selectedChannels, 'String')), 2)
    errordlg('Please select two input channels from ''Available Channels''.','Setting Error','modal')
    return;
end

%  Process Sanity check ( only check underlying data )
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
try
    userData.crtProc.sanityCheck;
catch ME
    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Retrieve GUI-defined parameters
channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
