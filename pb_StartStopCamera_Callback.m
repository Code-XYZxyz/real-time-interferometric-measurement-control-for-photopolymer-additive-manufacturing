% --- Executes on button press in pb_StartStopCamera.
function pb_StartStopCamera_Callback(hObject, eventdata, handles)
% hObject    handle to pb_StartStopCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global g_acquiring
global g_previewFrameIdx

handles = updateParameters(handles);

%% Start/Stop Camera
if strcmp(get(handles.pb_StartStopCamera,'String'),'Start Camera')
    % Camera is off. Change button string and start camera.
    set(handles.pb_StartStopCamera,'String','Stop Camera');
    set(handles.pb_StartStopCamera,'Enable','on');
    set(handles.pb_SetROI,'Enable','on');
    set(handles.pb_AcquireAVI,'String','Acquire & Analyze');
    set(handles.pb_AcquireAVI,'Enable','off');
    set(handles.pb_Snapshot,'Enable','on');
    set(handles.pb_OpenCloseUV,'Enable','off');
    set(handles.st_InterferogramStatusBar,'String','Starting Camera');

    %% Create video object
    ...Create and configure a video input object for the acquisition
        ...Choose one of the following 3 adaptors
        %---- Adpator 1: adaptor on my laptop
%     handles.video = videoinput('winvideo',1,'YUY2_640x480');
%     handles.video = videoinput('winvideo',1,'I420_640x480');
    
    % %---- Adaptor 2: incorrect adaptor on the lab computer "winvideo"
    % % handles.video = videoinput('winvideo',1,'YUY2_648x486');
    % % handles.video = videoinput('winvideo',1,'YUY2_640x480');
    % handles.video = videoinput('winvideo',1,'YUY2_2592x1944');
    % % Configure vidObj1's video source properties.
    % handles.srcObj1 = get(handles.video, 'Source');
    % set(handles.srcObj1(1), 'Brightness', 45);
    % set(handles.srcObj1(1), 'Exposure', -13);
    % set(handles.srcObj1(1), 'Gain', 3);
    
    %---- Adaptor 3: CORRECT adaptor "Gige" on the lab computer 648x486
    handles.video = videoinput('gige', 1, 'Mono8');
    handles.src = getselectedsource(handles.video);
    handles.src.BinningHorizontal = 4;
    handles.src.BinningVertical = 4;
    handles.src.ExposureTimeAbs = 105;
    handles.cp.FPS = str2double(get(handles.ed_FPS,'String')); % Acuiqistion Frames per Second
    handles.src.AcquisitionFrameRateAbs = handles.cp.FPS;
    
    % %---- Adaptor 4: adaptor on my home desktop
    % handles.video = videoinput('winvideo',1,'I420_640x480');
    
    %% ========= end of configuring video source
    
    handles.video.FramesPerTrigger = 1; % Acquire only one frame each time
    % handles.video.FramesPerTrigger = Inf; % keep acquiring frames until a stop command
    handles.video.TriggerRepeat = Inf; % Go on forever until stopped
    handles.video.ReturnedColorspace = 'grayscale';
    triggerconfig(handles.video,'manual');
    
    % image resolution, width, height
    vidRes = get(handles.video, 'VideoResolution');
    handles.cp.imW = vidRes(1); % image width
    handles.cp.imH = vidRes(2); % image height
    handles.cp.nBands = get(handles.video, 'NumberOfBands'); % # of bands, 1 for grayscale, 3 for RGB
    handles.cp.maxNumCache = 100;
    
    % is realtime
    handles.cp.isRT = 1;

    % flag if ROI is set or not
    handles.mask = 0;
    handles.toggle_POI_center = 1;
    handles.toggle_POI_corner = 0;
    handles.toggle_POI_horizon_line = 0;
    handles.toggle_POI_vertical_line = 0;
    handles.cp.POI = [];

    g_previewFrameIdx = 0;
    g_acquiring = 0;

    % use the timer to process input frames
    handles.video.TimerPeriod = 1 / handles.cp.FPS;
    handles.video.TimerFcn = {@Real_Time_ICM_processMeasureTimer, hObject};
    
    % video data logging setting
    %-- logging method 1:
    handles.video.LoggingMode = 'memory';

%     %-- logging method 2:
%     handles.video.LoggingMode = 'disk&memory';
%     diskLogger = VideoWriter(fullfile(handles.cp.VideoFolder,sprintf('ICM_%s.avi',datestr(now,'yyyymmdd_HHMMSS'))), 'Grayscale AVI');
% %     diskLogger = VideoWriter(fullfile(handles.cp.VideoFolder,sprintf('ICM_%s.avi',datestr(now,'yyyymmdd_HHMMSS'))), 'Uncompressed AVI');
%     diskLogger.FrameRate = handles.cp.FPS;
%     handles.video.DiskLogger = diskLogger;

    start(handles.video);
    guidata(hObject, handles);
    set(handles.st_InterferogramStatusBar,'String','Click "Set ROI" to proceed');
else
    % Camera is on. Stop camera and change button string.
    if isrunning(handles.video)
        stop(handles.video);
        flushdata(handles.video);%clear memory buffer
    end

    delete(handles.video)
    clear handles.video

    set(handles.pb_StartStopCamera,'String','Start Camera');
    set(handles.pb_AcquireAVI,'String','Acquire & Analyze');
    set(handles.pb_AcquireAVI,'Enable','off');
    set(handles.pb_Snapshot,'Enable','off');
    set(handles.pb_OpenCloseUV,'Enable','on');

    guidata(hObject, handles);
end

