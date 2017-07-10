% --- Executes on button press in pb_AcquireAVI.
function pb_AcquireAVI_Callback(hObject, eventdata, handles)
% hObject    handle to pb_AcquireAVI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global g_all_frame
global g_uvStatus
global g_acquiring
global g_POI
global g_mmf
global g_frameIdx
global g_job

handles = updateParameters(handles);

if strcmp(get(handles.pb_AcquireAVI,'String'), 'Acquire & Analyze')
    % Start acquisition
    nPOI = size(g_POI, 2);
    if nPOI == 0
        disp('ERROR: POI is empty!')
        set(handles.st_InterferogramStatusBar,'String','ERROR: POI is empty! Click "Set ROI" first to proceed!');
        return
    end

    % Camera is not acquiring. Change button string and start acquisition.
    set(handles.pb_AcquireAVI,'String','Stop Acquisition');
    set(handles.st_InterferogramStatusBar,'String','Acquiring will begin soon after ICM worker is ready and UV opens ');
    set(handles.pb_StartStopCamera, 'Enable','off');
    set(handles.pb_SetROI,'Enable','off');
    pause(0.5);

    % dynamic params
    handles.cp.POI = g_POI;
    disp(handles.cp)
    
%     % start pool for 1 main worker + N fitter
%     handles.cp.numFitter = 2;
%     poolobj = gcp('nocreate');
%     if isempty(poolobj)
%         numPool = handles.numFitter+1;
%         parpool(numPool);
%     end

    % init memory map file for gui session
    set(handles.st_InterferogramStatusBar,'String','creating memory map file');
    pause(1);
    clearFile = 1;
    g_mmf = icm_init_mem_file(handles.cp, clearFile);    

    %% Start the icm main worker
    disp('Starting ICM main worker')
    set(handles.st_InterferogramStatusBar,'String','Starting ICM main worker');
    pause(0.5);
    g_job = batch(@icm_main_worker, 1, {handles.cp},...
        'Pool', handles.cp.numFitter,...
        'AttachedFiles', {'icmFit2.m'}...
        );
    wait(g_job, 'running');
    disp('ICM main worker started, wait until worker is ready')
    set(handles.st_InterferogramStatusBar,'String','ICM main worker started, wait until worker is ready');
    while g_mmf.Data(1).status(3) == 0
        pause(1)
    end
    disp('icm main worker is ready!')
    set(handles.st_InterferogramStatusBar,'String','ICM main worker is ready! Click "Open UV" to start ECPL and ICM&M.');


    g_frameIdx = 0;
    g_acquiring = 1;
    g_uvStatus = 0;
    g_all_frame = zeros(handles.cp.imH, handles.cp.imW, 1, 1000); % HOME computer & camera

    set(handles.pb_OpenCloseUV,'Enable','on');
    
else
    % Stop acquisition
    % Camera is acquiring. Stop acquisition, save video data,
    % and change button string.
    % stop acquisition by stopping "trigger" in the Real_Time_ICM_processMeasureTimer.m"
%     set(handles.pb_AcquireAVI,'String','Acquire & Analyze');
    set(handles.pb_AcquireAVI,'Enable','off');
    set(handles.pb_StartStopCamera,'Enable','on');
    set(handles.pb_OpenCloseUV,'Enable','off');

    % if real-time, stop acquisition and record the video data to disk
    if strcmp(g_job.State, 'running') == 1
        % tell the worker the measurement stops
        disp('mark mmf to stop measurement by user action');
        set(handles.st_InterferogramStatusBar,'String','mark mmf to stop measurement by user action');
        g_mmf.Data(1).status(1) = 0;
    end

    g_acquiring = 0;

    %% write video file
    set(handles.st_InterferogramStatusBar,'String','Writing and saving ICM&M video file');
    nFrame = g_frameIdx;
    % Construct a VideoWriter object, which creates a Motion-JPEG AVI file by default.
    videoFile = fullfile(handles.cp.VideoFolder,sprintf('ICM_%s.avi',datestr(now,'yyyymmdd_HHMMSS')));
    fprintf('Start writing video file %s, nFrame %d\n', videoFile, nFrame);
    outputVideo = VideoWriter(videoFile, 'Grayscale AVI');
    outputVideo.FrameRate = handles.cp.FPS;
    open(outputVideo);
    %Loop through the image sequence, load each image, and then write it to the video.
    writeVideo(outputVideo, uint8(g_all_frame(:,:,:,1:nFrame))); % HOME computer & camera
    % Finalize the video file.
    close(outputVideo);
    fprintf('Wrote %d frames to video file %s\n', nFrame, videoFile);
    
%     % save metadata
%     S.frameIdx  = g_frameIdx;
%     S.cp = handles.cp;
%     save(strcat(handles.cp.ResultFolder,strcat('\IMAQ_',datestr(now,'yyyymmdd_HHMMSS'),'.mat')), '-struct', 'S');

    %% wait main worker done
    disp('wait for main worker done.');
    set(handles.st_InterferogramStatusBar,'String','Acquisition is done. Waiting for main worker to finish measurement and analysis.');
    wait(g_job)
    diary(g_job)
    handles.jobret = fetchOutputs(g_job);
    disp('returned result from main worker');
    set(handles.st_InterferogramStatusBar,'String','Main worker Returned measurement and analysis results.');
    disp(handles.jobret);
%     delete(g_job)
    assignin('base', 'ret', handles.jobret);
    g_mmf = 1;
    clear g_mmf;
    
end
guidata(hObject, handles);

