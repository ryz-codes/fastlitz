% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function varargout = fastlitz(varargin)
% FASTLITZ MATLAB code for fastlitz.fig
%      FASTLITZ, by itself, creates a new FASTLITZ or raises the existing
%      singleton*.
%
%      H = FASTLITZ returns the handle to a new FASTLITZ or the handle to
%      the existing singleton*.
%
%      FASTLITZ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FASTLITZ.M with the given input arguments.
%
%      FASTLITZ('Property','Value',...) creates a new FASTLITZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fastlitz_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fastlitz_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fastlitz

% Last Modified by GUIDE v2.5 26-Nov-2013 17:39:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fastlitz_OpeningFcn, ...
                   'gui_OutputFcn',  @fastlitz_OutputFcn, ...
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


% --- Executes just before fastlitz is made visible.
function fastlitz_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fastlitz (see VARARGIN)

% Choose default command line output for fastlitz
handles.output = hObject;

% Update display
handles = disp_refresh(handles);
handles = update_discret(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fastlitz wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fastlitz_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function txtDiam_Callback(hObject, eventdata, handles)
% hObject    handle to txtDiam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = disp_refresh(handles);
% Save data
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of txtDiam as text
%        str2double(get(hObject,'String')) returns contents of txtDiam as a double


% --- Executes during object creation, after setting all properties.
function txtDiam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtDiam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtDiamInsu_Callback(hObject, eventdata, handles)
% hObject    handle to txtDiamInsu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = disp_refresh(handles);
% Save data
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of txtDiamInsu as text
%        str2double(get(hObject,'String')) returns contents of txtDiamInsu as a double


% --- Executes during object creation, after setting all properties.
function txtDiamInsu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtDiamInsu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtStrands_Callback(hObject, eventdata, handles)
% hObject    handle to txtStrands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = disp_refresh(handles);
% Save data
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of txtStrands as text
%        str2double(get(hObject,'String')) returns contents of txtStrands as a double


% --- Executes during object creation, after setting all properties.
function txtStrands_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtStrands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtPitch_Callback(hObject, eventdata, handles)
% hObject    handle to txtPitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = disp_refresh(handles);
% Save data
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of txtPitch as text
%        str2double(get(hObject,'String')) returns contents of txtPitch as a double


% --- Executes during object creation, after setting all properties.
function txtPitch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtLen_Callback(hObject, eventdata, handles)
% hObject    handle to txtLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = disp_refresh(handles);
% Save data
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of txtLen as text
%        str2double(get(hObject,'String')) returns contents of txtLen as a double


% --- Executes during object creation, after setting all properties.
function txtLen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Collects the data from the text boxes, performs cleanups and error checks
function h = disp_refresh(h)
h.wlen = str2double(get(h.txtLen,'String'));
h.srad = str2double(get(h.txtDiam,'String'))/2;
insu = str2double(get(h.txtDiamInsu,'String'))/2;
h.insu = insu/h.srad - 1;
h.pitch = str2num(get(h.txtPitch,'String'));
h.ns = str2num(get(h.txtStrands,'String'));

% Generate visual model
% Detect if we want to draw a segment or a full
drawfull = get(h.optDraw,'value');
if drawfull
    seglen = h.wlen;
    set(h.optDraw,'String','Draw Full');
else
    seglen = 3*h.srad*3^numel(h.ns);
    set(h.optDraw,'String','Draw Seg');
end
dmod = litz.wire(seglen,h.pitch,h.ns,h.srad,h.insu);

cla(h.axes1);

% Show model
litz.showGeom(dmod,h.axes1);


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% siz = get(hObject,'position');
% sizax = get(handles.axes1,'position');
% sizax(3) = siz(3) - sizax(1);
% sizax(4) = siz(4) - sizax(2);
% set(handles.axes1,'position',sizax);


% --- Executes on button press in optDraw.
function optDraw_Callback(hObject, eventdata, handles)
% hObject    handle to optDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = disp_refresh(handles);
% Save data
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of optDraw


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = handles;
dmod = litz.wire(h.wlen,h.pitch,h.ns,h.srad,h.insu,h.dval,h.odd);
model = pmag(dmod);
assignin('base','model',model)
fprintf('Object ''model'' added to workspace.\nCommand line: \n\t[Z] = model.impedance(f) for impedance, or \n\t[P] = model.power(f,[Bx,By,Bz]) for B excitation, or\n\t[P] = model.power(f,[Hx,Hy,Hz],''A/m'') for H excitation\n');

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function txtDiscret_Callback(hObject, eventdata, handles)
% hObject    handle to txtDiscret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = update_discret(handles);
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of txtDiscret as text
%        str2double(get(hObject,'String')) returns contents of txtDiscret as a double


% --- Executes during object creation, after setting all properties.
function txtDiscret_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtDiscret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = update_discret(handles);
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox1

function h = update_discret(h)
h.dval = str2double(get(h.txtDiscret,'String'));
h.odd = get(h.checkbox1,'Value');
path = [0 0 0; 1 0 0];
wrad = 1;
mod = litz.solidify(path,wrad,h.dval,h.odd);
litz.showPEEC(mod,h.axes2);
view(h.axes2,120,20);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Copyright (c) 2013, Richard Y. Zhang
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its 
% contributors may be used to endorse or promote products derived from 
% this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
