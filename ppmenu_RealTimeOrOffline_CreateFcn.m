%% --- Executes during object creation, after setting all properties.
function ppmenu_RealTimeOrOffline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppmenu_RealTimeOrOffline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
