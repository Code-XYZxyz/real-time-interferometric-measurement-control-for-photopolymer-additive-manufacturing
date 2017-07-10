% --- Executes on selection change in ppmenu_RealTimeOrOffline.
function ppmenu_RealTimeOrOffline_Callback(hObject, eventdata, handles)
% hObject    handle to ppmenu_RealTimeOrOffline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppmenu_RealTimeOrOffline contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppmenu_RealTimeOrOffline
global ICMvidData
contents = get(hObject, 'Value');
switch contents
    case 1 % Real Time Measurement
        % clear memory for new data
        ICMvidData = []; % clear previously loaded ICMvidData in workspace
      
        set(handles.pb_StartStopCamera,'Enable','on');
        set(handles.pb_Snapshot,'Enable','on');
        set(handles.pb_AcquireAVI,'Enable','off');
        set(handles.pb_SetROI,'Enable','off');
        
        set(handles.pb_LoadICMVideo, 'Enable', 'off');
        set(handles.tb_PlayStopICMVideo, 'Enable', 'off');
        set(handles.pb_ICMVideo1stFrame, 'Enable', 'off');
        set(handles.pb_ICMVideoPrevFrame, 'Enable', 'off');
        set(handles.pb_ICMVideoNextFrame, 'Enable', 'off');
        set(handles.pb_ICMVideoEndFrame, 'Enable', 'off');
        set(handles.st_DisplayedFrame,'String','Frame 0 of 0');
        set(handles.st_InterferogramStatusBar,'String','Start Camera to Acquire ICM Video for Real Time Measurement');
    case 2 % Offline Calibration
        set(handles.pb_StartStopCamera,'Enable','off');
        set(handles.pb_Snapshot,'Enable','off');
        set(handles.pb_AcquireAVI,'Enable','off');
        set(handles.pb_SetROI,'Enable','on');
        
        set(handles.pb_LoadICMVideo, 'Enable', 'on');
        set(handles.tb_PlayStopICMVideo, 'Enable', 'on');
        set(handles.pb_ICMVideo1stFrame, 'Enable', 'on');
        set(handles.pb_ICMVideoPrevFrame, 'Enable', 'on');
        set(handles.pb_ICMVideoNextFrame, 'Enable', 'on');
        set(handles.pb_ICMVideoEndFrame, 'Enable', 'on');
        set(handles.st_InterferogramStatusBar,'String','Load an ICM Video for Offline Calibration');
    otherwise
end
guidata(hObject, handles);
