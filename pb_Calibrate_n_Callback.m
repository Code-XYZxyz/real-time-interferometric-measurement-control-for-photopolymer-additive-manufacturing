% --- Executes on button press in pb_Calibrate_n.
function pb_Calibrate_n_Callback(hObject, eventdata, handles)
%% ----- This function could:
...1. Calibrate refractive index using microscope measurement for future ICM measurement
...2. Save the calibrated refractive index n_m and delta_n


%% Calibrate
% hObject    handle to pb_Calibrate_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ICM_MeasuredHeight = str2double(get(handles.ed_ICM_MeasuredHeight,'String'));
handles.Phase2Pi = str2double(get(handles.ed_Phase2Pi,'String'));
handles.MicroscopeMeasHeight = str2double(get(handles.ed_MicroscopeMeasHeight,'String'));

% parameters
handles.FPS = str2double(get(handles.ed_FPS,'String'));
handles.MeasPeriodSamples = str2double(get(handles.ed_MeasPeriodSamples,'String'));
handles.Wavelength = str2double(get(handles.ed_Wavelength,'String'));
handles.n_L = str2double(get(handles.ed_n_L,'String'));

% Delta_n = n_m - n_L = lamda*Ti*Sum(Freq)/2/MicroscopeMeasHeight, 
...where Ti=RollFreq/FPS=MeasPeriodSample/FPS
handles.CalibratedDelta_n = handles.Wavelength*handles.Phase2Pi/2/handles.MicroscopeMeasHeight;
handles.Calibrated_n_m = handles.n_L + handles.CalibratedDelta_n;

set(handles.ed_CalibratedDelta_n,'string',handles.CalibratedDelta_n);
set(handles.ed_Calibrated_n_m,'string',handles.Calibrated_n_m);


guidata(hObject, handles);

%% Save calibration results
Calibration = struct('ICM_MeasuredHeight',handles.ICM_MeasuredHeight,'Phase2Pi',handles.Phase2Pi,...
    'MicroscopeMeasHeight',handles.MicroscopeMeasHeight,...
    'CalibratedDelta_n',handles.CalibratedDelta_n,'Calibrated_n_m',handles.Calibrated_n_m,...
    'Wavelength',handles.Wavelength,'n_L',handles.n_L,...
    'FPS',handles.FPS,'MeasPeriodSamples',handles.MeasPeriodSamples);
    
save(strcat(handles.cp.ResultFolder,strcat('\Calibration_',datestr(now,'yyyymmdd_HHMMSS'),'.mat')),'Calibration');
