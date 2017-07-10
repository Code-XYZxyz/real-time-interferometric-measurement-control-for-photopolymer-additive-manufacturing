% --- Executes when user attempts to close Real_Time_ICM.
function Real_Time_ICM_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Real_Time_ICM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

delete(imaqfind); % xyz inserted code

clear g_mmf;

%% ------------ delete parallel pool---------
p = gcp('nocreate');
if ~isempty(p)
    delete(p)
end;

%% ------------end of xyz inserted codes---------------
