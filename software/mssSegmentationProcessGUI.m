function varargout = mssSegmentationProcessGUI(varargin)
% MSSSEGMENTATIONPROCESSGUI M-file for mssSegmentationProcessGUI.fig
%      MSSSEGMENTATIONPROCESSGUI, by itself, creates a new MSSSEGMENTATIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = MSSSEGMENTATIONPROCESSGUI returns the handle to a new MSSSEGMENTATIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      MSSSEGMENTATIONPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSSSEGMENTATIONPROCESSGUI.M with the given input arguments.
%
%      MSSSEGMENTATIONPROCESSGUI('Property','Value',...) creates a new MSSSEGMENTATIONPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mssSegmentationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mssSegmentationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

% Edit the above text to modify the response to help mssSegmentationProcessGUI

% Last Modified by GUIDE v2.5 28-Sep-2011 15:05:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mssSegmentationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mssSegmentationProcessGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before mssSegmentationProcessGUI is made visible.
function mssSegmentationProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)


processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Parameter setup
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
funParams = userData.crtProc.funParams_;

set(handles.listbox_Scales,'String',funParams.Scales);
value = find(funParams.FilterOrder==[1 3 5]);
set(handles.popupmenu_FilterOrder,'Value',value);

% Update user data and GUI data
handles.output = hObject;
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = mssSegmentationProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)


%  Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
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
funParams.Scales=str2num(get(handles.listbox_Scales,'String'));
props = get(handles.popupmenu_FilterOrder,{'String','Value'});
funParams.FilterOrder = str2double(props{1}{props{2}});

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

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
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on button press in pushbutton_addScale.
function pushbutton_addScale_Callback(hObject, eventdata, handles)

value = str2double(get(handles.edit_Scale,'String'));
if isnan(value) || value<0
    errordlg('Please enter a valid scale for filtering');
    return;
else
    scales = str2num(get(handles.listbox_Scales,'String')); %#ok<ST2NM>
    if ismember(value,scales), return; end
    scales = sort([scales; value]);
    set(handles.listbox_Scales,'String',scales);
end

% --- Executes on button press in pushbutton_removeScale.
function pushbutton_removeScale_Callback(hObject, eventdata, handles)

props = get(handles.listbox_Scales,{'String','Value'});
if isempty(props{1}), return; end

props{1}(props{2})=[];
set(handles.listbox_Scales,'String',props{1},'Value',max(1,props{2}-1));
