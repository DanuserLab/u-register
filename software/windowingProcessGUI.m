function varargout = windowingProcessGUI(varargin)
% windowingProcessGUI M-file for windowingProcessGUI.fig
%      windowingProcessGUI, by itself, creates a new windowingProcessGUI or raises the existing
%      singleton*.
%
%      H = windowingProcessGUI returns the handle to a new windowingProcessGUI or the handle to
%      the existing singleton*.
%
%      windowingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in windowingProcessGUI.M with the given input arguments.
%
%      windowingProcessGUI('Property','Value',...) creates a new windowingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before windowingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to windowingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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

% Edit the above text to modify the response to help windowingProcessGUI

% Last Modified by GUIDE v2.5 01-Mar-2012 20:31:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @windowingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @windowingProcessGUI_OutputFcn, ...
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


% --- Executes just before windowingProcessGUI is made visible.
function windowingProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;
userData.numParams ={'ParaSize','PerpSize','MinSize','StartContour'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),...
    userData.numParams)
editSize(hObject,eventdata,handles);

% Read available segmentaation processes
segProc =  cellfun(@(x) isa(x,'MaskProcess'),userData.MD.processes_);
segProcID=find(segProc);
segProcNames = cellfun(@(x) x.getName(),userData.MD.processes_(segProc),'Unif',false);
segProcString = vertcat('Choose later',segProcNames(:));
segProcData=horzcat({[]},num2cell(segProcID));

% Read the default segmentation process index
% If empty, try to propagate segmentation process from protrusion process

ProtProcIndex = find(cellfun(@(x) isa(x,'ProtrusionProcess'),userData.crtPackage.processes_));

initSegProcIndex = funParams.SegProcessIndex;
if isempty(initSegProcIndex) && ~isempty(ProtProcIndex)
    if ~isempty(userData.crtPackage.processes_{ProtProcIndex}.funParams_.SegProcessIndex)
        initSegProcIndex = userData.crtPackage.processes_{ProtProcIndex}.funParams_.SegProcessIndex;  
    end
end
segProcValue = find(cellfun(@(x) isequal(x,initSegProcIndex),segProcData));
if isempty(segProcValue), segProcValue = 1; end
set(handles.popupmenu_SegProcessIndex,'String',segProcString,...
    'UserData',segProcData,'Value',segProcValue);

% Create pop-up menu for windowing methods
methods = WindowingProcess.getMethods;
methodValue = find(strcmpi(funParams.MethodName,{methods.name}));
set(handles.popupmenu_MethodName,'String',{methods.description},...
    'UserData',{methods.name},'Value',methodValue);

% Update reinitialization parameters
if isinf(funParams.ReInit)
    set(handles.checkbox_doReInit,'Value',0);
    set(handles.edit_ReInit,'String','','Enable','off');
else
    set(handles.checkbox_doReInit,'Value',1);
    set(handles.edit_ReInit,'String',funParams.ReInit,'Enable','on');
end

% Update StartPoint
userData.StartPoint=funParams.StartPoint;

set(handles.checkbox_useStartPoint,'Value',~isempty(funParams.StartPoint));
if isempty(funParams.StartPoint)
    set(handles.edit_StartPoint,'String',funParams.StartPoint);
else
    set(handles.edit_StartPoint,'String',['(' num2str(funParams.StartPoint(1)) ',' num2str(funParams.StartPoint(2)) ')']);
end

userData.previewFig=-1;
userData.imPointHandle.isvalid=0;

% Update channels listboxes depending on the selected process
popupmenu_SegProcessIndex_Callback(hObject, eventdata, handles)

% Choose default command line output for windowingProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = windowingProcessGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

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


function editSize(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
pixelSize = userData.MD.pixelSize_;
paraSize = str2double(get(handles.edit_ParaSize,'String'));
perpSize = str2double(get(handles.edit_PerpSize,'String'));
paraSizeMicrons = pixelSize/1000*paraSize;
perpSizeMicrons = pixelSize/1000*perpSize;
if ~isempty(paraSizeMicrons) && ~isnan(paraSizeMicrons)
    set(handles.edit_ParaSizeMicrons,'String',paraSizeMicrons);
else
    set(handles.edit_ParaSizeMicrons,'String','');
end
if ~isempty(perpSizeMicrons) && ~isnan(perpSizeMicrons)
    set(handles.edit_PerpSizeMicrons,'String',perpSizeMicrons);
else
    set(handles.edit_PerpSizeMicrons,'String','');
end


% --- Executes on button press in checkbox_doReInit.
function checkbox_doReInit_Callback(hObject, eventdata, handles)

if ~get(hObject,'Value')
    set(handles.edit_ReInit,'String','','Enable','off');
else
    set(handles.edit_ReInit,'Enable','on');
end

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Check user input
userData = get(handles.figure1, 'UserData');
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
else
    channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
    funParams.ChannelIndex = channelIndex;
end

for i=1:numel(userData.numParams)
    value = str2double(get(handles.(['edit_' userData.numParams{i}]),'String'));
    if ~isposint(value)
        errordlg(['Please enter a valid value for the '...
            get(handles.(['text_' userData.numParams{i}]),'String') '.'],'Setting Error','modal');
        return;
    end
    funParams.(userData.numParams{i}) = value; 
end

% Retrieve mask process index and class (for propagation)
props=get(handles.popupmenu_SegProcessIndex,{'UserData','Value'});
funParams.SegProcessIndex = props{1}{props{2}};
if ~isempty(funParams.SegProcessIndex)
    segProcessClass=class(userData.MD.processes_{funParams.SegProcessIndex});
else
    segProcessClass = '';
end
% Retrieve windowing method
props=get(handles.popupmenu_MethodName,{'UserData','Value'});
funParams.MethodName=props{1}{props{2}};

if get(handles.checkbox_doReInit,'Value')
    reInit = str2double(get(handles.edit_ReInit,'String'));
    if isnan(reInit) || reInit<1
        errordlg(['Please enter a valid value for the '...
            get(handles.text_ReInit,'String') '.'],'Setting Error','modal');
        return;
    end
    funParams.ReInit=reInit; 
else
    funParams.ReInit=Inf;
end

if get(handles.checkbox_useStartPoint,'Value')
    if userData.imPointHandle.isvalid
        userData.StartPoint=ceil(getPosition(userData.imPointHandle));
    end
    funParams.StartPoint=userData.StartPoint;
    if ishandle(userData.previewFig), delete(userData.previewFig); end
end
% Process Sanity check ( only check underlying data )
try
    userData.crtProc.sanityCheck;
catch ME
    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters
setMaskProcess = @(x) parseProcessParams(x, struct('SegProcessIndex',...
    x.owner_.getProcessIndex(segProcessClass,1,false)));
processGUI_ApplyFcn(hObject, eventdata, handles,funParams,'settingFcn',{setMaskProcess});


% --- Executes on button press in checkbox_selectStartPoint.
function update_data(hObject, eventdata, handles)
userData=get(handles.figure1,'UserData');
if ~get(handles.checkbox_selectStartPoint,'Value')
    if userData.imPointHandle.isvalid,
        userData.StartPoint=ceil(getPosition(userData.imPointHandle));
        userData.StartPoint=userData.StartPoint;
    end
    if ishandle(userData.previewFig), delete(userData.previewFig); end
else
    if ishandle(userData.previewFig), 
        delete(userData.previewFig); 
    else
        userData.previewFig=figure('NumberTitle','off','Name','Select the origin of the windows',...
            'DeleteFcn',@close_previewFig,'UserData',handles.figure1);
    end
    
    imHandle = findobj(userData.previewFig,'Type','image');
    if isempty(imHandle)
        selectedChannels = get(handles.listbox_selectedChannels,'UserData');
        selectedChannels(4:end)=[];
        imHandle=userData.MD.channels_(selectedChannels).draw(1);
    end

    
    if userData.imPointHandle.isvalid
        setPosition(userData.imPointHandle,userData.StartPoint);
    else
        % Create a new imrect object and store the handle
        if isempty(userData.StartPoint),
            userData.StartPoint=userData.MD.imSize_(2:-1:1)/2;
        end
        userData.imPointHandle = impoint(get(imHandle,'Parent'),userData.StartPoint);
        constraintFcn = makeConstrainToRectFcn('impoint',get(imHandle,'XData'),get(imHandle,'YData'));
        motionFcn = @(x) updateStartPointPosition(x,userData.previewFig);
        setPositionConstraintFcn(userData.imPointHandle,constraintFcn);
        addNewPositionCallback(userData.imPointHandle,motionFcn);
        updateStartPointPosition(userData.StartPoint,userData.previewFig)
    end
    
end
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

function close_previewFig(hObject, eventdata)
handles = guidata(get(hObject,'UserData'));
set(handles.checkbox_selectStartPoint,'Value',0);
update_data(handles.checkbox_selectStartPoint, eventdata, handles);

function updateStartPointPosition(pos,fig)
% Pos is in xy coordinate while start point is in
handles = guidata(get(fig,'UserData'));
set(handles.edit_StartPoint,'String',['(' num2str(ceil(pos(1))) ',' num2str(ceil(pos(2))) ')']);
       
