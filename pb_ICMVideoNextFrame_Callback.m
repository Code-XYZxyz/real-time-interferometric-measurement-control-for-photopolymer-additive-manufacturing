% --- Executes on button press in pb_ICMVideoNextFrame.
function pb_ICMVideoNextFrame_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ICMVideoNextFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ICMvidData
global hImage

if isempty(ICMvidData),return, end;
if handles.CFrameInd ~= handles.AVInFrame % not end frame currently
    handles.CFrameInd = handles.CFrameInd + 1; % forward to next frame
end
CFrame=ICMvidData(:,:,handles.CFrameInd);
% axes(handles.Interferogram);
% imshow(CFrame);
hImage = imshow(CFrame, 'Parent', handles.Interferogram);
set(handles.st_DisplayedFrame,'String',sprintf('Frame %d of %d',handles.CFrameInd,handles.AVInFrame));
guidata(hObject, handles);  
