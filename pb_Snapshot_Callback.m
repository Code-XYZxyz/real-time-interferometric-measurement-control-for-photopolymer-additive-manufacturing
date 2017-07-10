% --- Executes on button press in pb_Snapshot.
function pb_Snapshot_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Snapshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% ----------beginning of xyz inserted code----------
if ~isrunning(handles.video), return,end; % if camera not started yet
% Capturing a single image and Save to current location
snappedFrame = getsnapshot(handles.video);
FileNameSnapFrame = strcat('Interferogram_',datestr(now,'yyyymmdd_HHMMSS'));
imwrite (snappedFrame,fullfile(handles.cp.VideoFolder,strcat(FileNameSnapFrame,'.jpg')),'jpg');
save(fullfile(handles.cp.VideoFolder,strcat(FileNameSnapFrame,'.mat')),'snappedFrame');
set(handles.st_InterferogramStatusBar,'String',sprintf('Snapshot saved to: %s',strcat(FileNameSnapFrame,'.mat and .jpg')));

% % Capture and Save the current displayed frame
% CurrentDisplayedFrame = get(get(handles.Interferogram,'children'),'cdata'); 
% FileNameCurrentFrame = strcat('Interferogram_',datestr(now,'yyyymmdd_HHMMSS'));
% imwrite (CurrentDisplayedFrame,strcat(FileNameCurrentFrame,'.jpg'),'jpg');
% save(strcat(FileNameCurrentFrame,'.mat'), 'CurrentDisplayedFrame');
% disp(strcat('Current Displayed Frame saved to file:  ',strcat(FileNameCurrentFrame,'.mat')));
% ------------end of xyz inserted codes---------------
