function varargout = TremblingCorrectionProcessGUI(varargin)
%TREMBLINGCORRECTIONPROCESSGUI MATLAB code file for TremblingCorrectionProcessGUI.fig
%      TREMBLINGCORRECTIONPROCESSGUI, by itself, creates a new TREMBLINGCORRECTIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = TREMBLINGCORRECTIONPROCESSGUI returns the handle to a new TREMBLINGCORRECTIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      TREMBLINGCORRECTIONPROCESSGUI('Property','Value',...) creates a new TREMBLINGCORRECTIONPROCESSGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to TremblingCorrectionProcessGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TREMBLINGCORRECTIONPROCESSGUI('CALLBACK') and TREMBLINGCORRECTIONPROCESSGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TREMBLINGCORRECTIONPROCESSGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

% Edit the above text to modify the response to help TremblingCorrectionProcessGUI

% Last Modified by GUIDE v2.5 12-Aug-2021 15:52:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TremblingCorrectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @TremblingCorrectionProcessGUI_OutputFcn, ...
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


% --- Executes just before TremblingCorrectionProcessGUI is made visible.
function TremblingCorrectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Parameter setup
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
funParams = userData.crtProc.funParams_;


% set GUI with Parameters
segProc =  cellfun(@(x) isa(x,'MaskProcess'),userData.MD.processes_);
segProcID=find(segProc);
segProcNames = cellfun(@(x) x.getName(),userData.MD.processes_(segProc),'Unif',false);
segProcString = vertcat('Choose later',segProcNames(:));
segProcData=horzcat({[]},num2cell(segProcID));
segProcValue = find(cellfun(@(x) isequal(x,funParams.SegProcessIndex),segProcData));
if isempty(segProcValue), segProcValue = 1; end
set(handles.popupmenu_SegProcessIndex,'String',segProcString,...
    'UserData',segProcData,'Value',segProcValue);

set(handles.edit_refineClosureRadius, 'String',num2str(funParams.refinementClosureRadius))

% Update channels listboxes depending on the selected process
popupmenu_SegProcessIndex_Callback(hObject, eventdata, handles)

% Update user data and GUI data
handles.output = hObject;
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = TremblingCorrectionProcessGUI_OutputFcn(hObject, eventdata, handles)
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

% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

% -------- Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

if isnan(str2double(get(handles.edit_refineClosureRadius, 'String'))) ...
    || str2double(get(handles.edit_refineClosureRadius, 'String')) < 0
  errordlg('Please provide a valid input for ''Refinement Closure Radius''.','Setting Error','modal');
  return;
end

% -------- Process Sanity check --------
% ( only check underlying data )
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

funParams.refinementClosureRadius = str2double(get(handles.edit_refineClosureRadius, 'String'));

% Retrieve mask process index and class (for propagation)
props=get(handles.popupmenu_SegProcessIndex,{'UserData','Value'});
funParams.SegProcessIndex = props{1}{props{2}};
if ~isempty(funParams.SegProcessIndex)
    segProcessClass=class(userData.MD.processes_{funParams.SegProcessIndex});
else
    segProcessClass = '';
end

% Set parameters and update main window
setMaskProcess = @(x) parseProcessParams(x, struct('SegProcessIndex',...
    x.owner_.getProcessIndex(segProcessClass,1,false)));
processGUI_ApplyFcn(hObject, eventdata, handles,funParams,'settingFcn',{setMaskProcess});


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end



% --- Executes on selection change in popupmenu_SegProcessIndex.
function popupmenu_SegProcessIndex_Callback(hObject, eventdata, handles)

% Retrieve selected process ID
props= get(handles.popupmenu_SegProcessIndex,{'UserData','Value'});
procID = props{1}{props{2}};

% Read process and check available channels
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
if isempty(procID)
    allChannelIndex=1:numel(userData.MD.channels_);
else
    allChannelIndex = find(userData.MD.processes_{procID}.checkChannelOutput);
end

% Set up available channels listbox
if ~isempty(allChannelIndex)
    if isempty(procID)
        channelString = userData.MD.getChannelPaths(allChannelIndex);
    else
        channelString = userData.MD.processes_{procID}.outFilePaths_(1,allChannelIndex);
    end
else
    channelString = {};
end
set(handles.listbox_availableChannels,'String',channelString,'UserData',allChannelIndex);

% Set up selected channels listbox
channelIndex = get(handles.listbox_selectedChannels, 'UserData');
channelIndex = intersect(channelIndex,allChannelIndex);
if ~isempty(channelIndex)
    if isempty(procID)
        channelString = userData.MD.getChannelPaths(channelIndex);
    else
        channelString = userData.MD.processes_{procID}.outFilePaths_(1,channelIndex);
    end
else
    channelString = {};
end
set(handles.listbox_selectedChannels,'String',channelString,'UserData',channelIndex);
