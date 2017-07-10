% --- Executes on button press in pb_ICMVideo1stFrame.
function pb_ICMVideo1stFrame_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ICMVideo1stFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ICMvidData
global hImage

if isempty(ICMvidData),return, end;
handles.CFrameInd = 1; % show first frame
CFrame=ICMvidData(:,:,handles.CFrameInd);
% axes(handles.Interferogram);
% imshow(CFrame);
hImage = imshow(CFrame, 'Parent', handles.Interferogram);
set(handles.st_DisplayedFrame,'String',sprintf('Frame %d of %d',handles.CFrameInd,handles.AVInFrame));
guidata(hObject, handles);
