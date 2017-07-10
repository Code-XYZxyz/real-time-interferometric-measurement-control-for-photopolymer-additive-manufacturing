function varargout = Real_Time_ICM(varargin)
% REAL_TIME_ICM MATLAB code for Real_Time_ICM.fig
%      REAL_TIME_ICM, by itself, creates a new REAL_TIME_ICM or raises the existing
%      singleton*.
%
%      H = REAL_TIME_ICM returns the handle to a new REAL_TIME_ICM or the handle to
%      the existing singleton*.
%
%      REAL_TIME_ICM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REAL_TIME_ICM.M with the given input arguments.
%
%      REAL_TIME_ICM('Property','Value',...) creates a new REAL_TIME_ICM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Real_Time_ICM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Real_Time_ICM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Real_Time_ICM

% Last Modified by GUIDE v2.5 25-Jul-2016 12:01:00

%% clear persistent values
clear Real_Time_ICM_processMeasureTimer;
clear icm_main_worker;
% clear global;
% clear -regexp ^g_;

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Real_Time_ICM_OpeningFcn, ...
                   'gui_OutputFcn',  @Real_Time_ICM_OutputFcn, ...
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

function ed_FPS_Callback(hObject, eventdata, handles)
handles.FPS = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_FPS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_MeasPeriodSamples_Callback(hObject, eventdata, handles)
handles.MeasPeriodSamples = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_MeasPeriodSamples_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_MovingHorizonL_Callback(hObject, eventdata, handles)
handles.MovingHorizonL = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_MovingHorizonL_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_HalfLife_Callback(hObject, eventdata, handles)
handles.HalfLife = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_HalfLife_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% function ed_GOF_rSqaure_Callback(hObject, eventdata, handles)
% handles.GOF_rSquare = str2double(get(hObject,'String'));
% guidata(hObject, handles);
% % --- Executes during object creation, after setting all properties.
% function ed_GOF_rSqaure_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

function ed_Wavelength_Callback(hObject, eventdata, handles)
handles.Wavelength = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_Wavelength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_n_L_Callback(hObject, eventdata, handles)
handles.n_L = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_n_L_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% 
function ed_n_m_Callback(hObject, eventdata, handles)
handles.n_m = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_n_m_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Begin:set the threhold to identify outlier frequency
 % 2 criterior used: freq<f_max && diff(freq)<f_diff_max 
% set the upper bound of the frequency change to remove the extrem outliers
function ed_f_diff_max_Callback(hObject, eventdata, handles)
handles.I_u = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_f_diff_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Pc: the threhold for identifying outliers in fitted coefficient of frequency
function ed_f_max_Callback(hObject, eventdata, handles)
handles.P_c = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_f_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% --- END:set the threhold to identify
...the beginning threhold period of zero-curing & curve fitting outliers

%% --- Begin:select the part of ROI to be measured: center, corner, centerline
% Executes on button press in checkbox_CenterPointROI.
function checkbox_CenterPointROI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CenterPointROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CenterPointROI
handles.toggle_POI_center = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes on button press in checkbox_CornerPointsROI.
function checkbox_CornerPointsROI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CornerPointsROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CornerPointsROI
handles.toggle_POI_corner = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes on button press in checkbox_hCenterlineROI.
function checkbox_hCenterlineROI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_hCenterlineROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_hCenterlineROI
handles.toggle_POI_horizon_line = get(hObject, 'Value');
guidata(hObject, handles);


% --- Executes on button press in checkbox_vCenterlineROI.
function checkbox_vCenterlineROI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vCenterlineROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vCenterlineROI
handles.toggle_POI_vertial_line = get(hObject, 'Value');
guidata(hObject, handles);
%% --- END:select the part of ROI to be measured: center, corner, centerline

function ed_GOF_rSquare_Callback(hObject, eventdata, handles)
handles.GOF_rSquare = str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ed_GOF_rSquare_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_GOF_rSquare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_ICM_MeasuredHeight_Callback(hObject, eventdata, handles)
% hObject    handle to ed_ICM_MeasuredHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_ICM_MeasuredHeight as text
%        str2double(get(hObject,'String')) returns contents of ed_ICM_MeasuredHeight as a double


% --- Executes during object creation, after setting all properties.
function ed_ICM_MeasuredHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_ICM_MeasuredHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_MicroscopeMeasHeight_Callback(hObject, eventdata, handles)
% hObject    handle to ed_MicroscopeMeasHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_MicroscopeMeasHeight as text
%        str2double(get(hObject,'String')) returns contents of ed_MicroscopeMeasHeight as a double


% --- Executes during object creation, after setting all properties.
function ed_MicroscopeMeasHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_MicroscopeMeasHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_Phase2Pi_Callback(hObject, eventdata, handles)
% hObject    handle to ed_Phase2Pi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_Phase2Pi as text
%        str2double(get(hObject,'String')) returns contents of ed_Phase2Pi as a double


% --- Executes during object creation, after setting all properties.
function ed_Phase2Pi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_Phase2Pi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_CalibratedDelta_n_Callback(hObject, eventdata, handles)
% hObject    handle to ed_CalibratedDelta_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_CalibratedDelta_n as text
%        str2double(get(hObject,'String')) returns contents of ed_CalibratedDelta_n as a double


% --- Executes during object creation, after setting all properties.
function ed_CalibratedDelta_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_CalibratedDelta_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_Calibrated_n_m_Callback(hObject, eventdata, handles)
% hObject    handle to ed_Calibrated_n_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_Calibrated_n_m as text
%        str2double(get(hObject,'String')) returns contents of ed_Calibrated_n_m as a double


% --- Executes during object creation, after setting all properties.
function ed_Calibrated_n_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_Calibrated_n_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in pb_Calibrate_n.
% function pb_Calibrate_n_Callback(hObject, eventdata, handles)
% % hObject    handle to pb_Calibrate_n (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)



function ed_SamplesNumBeforeMeasure_Callback(hObject, eventdata, handles)
% hObject    handle to ed_SamplesNumBeforeMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_SamplesNumBeforeMeasure as text
%        str2double(get(hObject,'String')) returns contents of ed_SamplesNumBeforeMeasure as a double


% --- Executes during object creation, after setting all properties.
function ed_SamplesNumBeforeMeasure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_SamplesNumBeforeMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_recWidth_Callback(hObject, eventdata, handles)
% hObject    handle to ed_recWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_recWidth as text
%        str2double(get(hObject,'String')) returns contents of ed_recWidth as a double


% --- Executes during object creation, after setting all properties.
function ed_recWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_recWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_recHeight_Callback(hObject, eventdata, handles)
% hObject    handle to ed_recHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_recHeight as text
%        str2double(get(hObject,'String')) returns contents of ed_recHeight as a double


% --- Executes during object creation, after setting all properties.
function ed_recHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_recHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_recGrayscale_Callback(hObject, eventdata, handles)
% hObject    handle to ed_recGrayscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_recGrayscale as text
%        str2double(get(hObject,'String')) returns contents of ed_recGrayscale as a double


% --- Executes during object creation, after setting all properties.
function ed_recGrayscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_recGrayscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_uvIris_Callback(hObject, eventdata, handles)
% hObject    handle to ed_uvIris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_uvIris as text
%        str2double(get(hObject,'String')) returns contents of ed_uvIris as a double


% --- Executes during object creation, after setting all properties.
function ed_uvIris_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_uvIris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_TargetCuredHeight_Callback(hObject, eventdata, handles)
% hObject    handle to ed_TargetCuredHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_TargetCuredHeight as text
%        str2double(get(hObject,'String')) returns contents of ed_TargetCuredHeight as a double


% --- Executes during object creation, after setting all properties.
function ed_TargetCuredHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_TargetCuredHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_ExposureTime_Callback(hObject, eventdata, handles)
% hObject    handle to ed_ExposureTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_ExposureTime as text
%        str2double(get(hObject,'String')) returns contents of ed_ExposureTime as a double


% --- Executes during object creation, after setting all properties.
function ed_ExposureTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_ExposureTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_ExposedCuredHeight_Callback(hObject, eventdata, handles)
% hObject    handle to ed_ExposedCuredHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_ExposedCuredHeight as text
%        str2double(get(hObject,'String')) returns contents of ed_ExposedCuredHeight as a double


% --- Executes during object creation, after setting all properties.
function ed_ExposedCuredHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_ExposedCuredHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_DarkCuredHeight_Callback(hObject, eventdata, handles)
% hObject    handle to ed_DarkCuredHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_DarkCuredHeight as text
%        str2double(get(hObject,'String')) returns contents of ed_DarkCuredHeight as a double


% --- Executes during object creation, after setting all properties.
function ed_DarkCuredHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_DarkCuredHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_TargetExpTime_Callback(hObject, eventdata, handles)
% hObject    handle to ed_TargetExpTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_TargetExpTime as text
%        str2double(get(hObject,'String')) returns contents of ed_TargetExpTime as a double


% --- Executes during object creation, after setting all properties.
function ed_TargetExpTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_TargetExpTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_NumFitters_Callback(hObject, eventdata, handles)
% hObject    handle to ed_NumFitters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_NumFitters as text
%        str2double(get(hObject,'String')) returns contents of ed_NumFitters as a double


% --- Executes during object creation, after setting all properties.
function ed_NumFitters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_NumFitters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
