% --- Executes on button press in pb_ICMVideoPrevFrame.
function pb_ICMVideoPrevFrame_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ICMVideoPrevFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ICMvidData
global hImage

if isempty(ICMvidData),return, end;
if handles.CFrameInd ~= 1 % not first frame currently
    handles.CFrameInd = handles.CFrameInd - 1; % back to previous frame
end
CFrame=ICMvidData(:,:,handles.CFrameInd);
% axes(handles.Interferogram);
% imshow(CFrame);
hImage = imshow(CFrame, 'Parent', handles.Interferogram);
set(handles.st_DisplayedFrame,'String',sprintf('Frame %d of %d',handles.CFrameInd,handles.AVInFrame));
guidata(hObject, handles);    
