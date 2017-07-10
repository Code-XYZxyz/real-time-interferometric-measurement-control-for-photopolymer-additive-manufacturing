%% process measurement timer func
function Real_Time_ICM_processMeasureTimer(vid, event, hObject)

global hImage
global g_all_frame
global uv
global g_uvStatus
global hdvp
global g_acquiring
global g_mmf
global g_frameIdx
global g_tImaqStart
global g_previewFrameIdx

%% Communication variables
global g_mmfWriteIdx

% persistent lastTic
% if ~isempty(lastTic)
%     timerDelay = toc(lastTic);
%     fprintf('delay=%.3f\n', timerDelay);
% end
% lastTic = tic;

%% Init everything and folders
handles = guidata(hObject);
if ~isrunning(vid)
    return
end

%% Real-time Preview video 
% if user does NOT push the button to acquire and log AVI, just preview
frame = getsnapshot(vid);
g_previewFrameIdx = g_previewFrameIdx + 1;
if g_previewFrameIdx == 1
    hImage = imshow(frame, 'Parent', handles.Interferogram);
else
    set(hImage, 'CData', frame);
end

if g_acquiring == 0 || g_uvStatus == 0
    return
end; 


if g_mmf.Data(1).status(1) == 0
    % capture stopped
    fprintf('Timer get msg that capture is stopped!\n')
    guidata(hObject, handles);
    return
end

if g_mmf.Data(1).status(3) == 2 && g_uvStatus == 1
    % shutdown uv
    icm_set_uv_status(0,handles.cp.uvIris);
    
    % display exposure time dynamically till UV closed
    set(handles.ed_ExposureTime, 'string', toc(g_tImaqStart)); 
    
    % Change button string
    set(handles.pb_OpenCloseUV, 'String', 'Open UV Light');
    set(handles.pb_OpenCloseUV, 'Enable', 'off');

    % update the Interferogram Status Bar for info
    set(handles.st_InterferogramStatusBar,'String','ON Target. Controller closed UV shutter. Click "Stop Acquisition" when ready.');

end


%% now we are acquiring and uvStatus is 1
if g_frameIdx == 0
%     g_tImaqStart = tic;
    fprintf('Measurement starts...\n')
    disp(handles.cp)
    g_mmfWriteIdx = 0;
end

%% Real-time acquisition:get latest frame and display
% if user has started to acquire AVI, start acquiring and logging data
% trigger(handles.video);
% frame = getdata(handles.video,1); % takes 27ms, slower than "getsnapshot"
% frame = getsnapshot(vid); % takes 7ms, much faster than "getdata"
% tImaqStart = tic; % start time of imaq
frameTime = toc(g_tImaqStart); % absolute elapsed time from start to current frame
g_frameIdx = g_frameIdx + 1;

%% process mmf
% read status

g_mmfWriteIdx = mod(g_mmfWriteIdx, handles.cp.maxNumCache) + 1;
% fprintf('g_mmfWriteIdx=%d\n', g_mmfWriteIdx)

% write frame
fprintf('Writing frame %d to cache %d\n', g_frameIdx, g_mmfWriteIdx)
if g_mmf.Data(g_mmfWriteIdx).frameIdx(1) ~= 0
    fprintf('Previous frame %d at cache %d is not processed yet!!!\n', g_mmf.Data(g_mmfWriteIdx).frameIdx(1), g_mmfWriteIdx)
end
g_mmf.Data(g_mmfWriteIdx).frame = frame;
g_mmf.Data(g_mmfWriteIdx).frameIdx(1) = g_frameIdx;
g_mmf.Data(g_mmfWriteIdx).frameTime(1) = frameTime;
g_mmf.Data(g_mmfWriteIdx).uvIris(1) = handles.cp.uvIris;
g_mmf.Data(g_mmfWriteIdx).snapTic(1) = tic;
g_mmf.Data(g_mmfWriteIdx).uvStatus(1) = g_uvStatus;

% update status
g_mmf.Data(1).status(2) = g_frameIdx;
set(hImage, 'CData', frame);
g_all_frame(:,:,1,g_frameIdx) = frame;

set(handles.st_DisplayedFrame,'String',sprintf('Frame %d',g_frameIdx));

% display exposure time dynamically till UV closed
if g_uvStatus == 1
    set(handles.ed_ExposureTime, 'string', frameTime); 
end

set(handles.ed_ICM_MeasuredHeight,'String',g_mmf.Data(1).avgTotalHeight(1));
set(handles.ed_Phase2Pi,'String',g_mmf.Data(1).avgTotalPhase(1));
%%
guidata(hObject, handles);

%% Real-time EC2C control 