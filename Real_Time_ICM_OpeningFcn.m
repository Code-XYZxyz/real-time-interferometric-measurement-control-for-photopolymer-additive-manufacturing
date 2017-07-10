% --- Executes just before Real_Time_ICM is made visible.
function Real_Time_ICM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Real_Time_ICM (see VARARGIN)

% Choose default command line output for Real_Time_ICM
handles.output = hObject;

% set the position of the GUI on screen (08-23-2016)
 set(handles.Real_Time_ICM,'Units', 'pixels');
 screenSize = get(0, 'ScreenSize');
 set(handles.Real_Time_ICM,'Position', [screenSize(1) screenSize(2) screenSize(3) screenSize(4)]);

%% Create a folder for results
rootDir = pwd;
handles.cp.ResultFolder = fullfile(rootDir, 'Result'); 
% create a result folder if it doesn't exist
if ~exist(handles.cp.ResultFolder, 'dir')
    mkdir(rootDir, 'Result');
end 

%% % Set up logging of acquired video file
% handles.video.LoggingMode = 'disk&memory';
% logging directory
handles.cp.VideoFolder = fullfile(handles.cp.ResultFolder, 'Video');
% create a video folder if it doesn't exist
if ~exist(handles.cp.VideoFolder, 'dir')
    mkdir(handles.cp.VideoFolder);
end
% % % handles.VideoLoggingNO = 1; % added field of video sequence NO.
% % diskLogger = VideoWriter(fullfile(handles.cp.VideoFolder,sprintf('ICM_%s.avi',datestr(now,'yyyymmdd_HHMMSS'))), 'Grayscale AVI');
% diskLogger = VideoWriter(fullfile(handles.cp.VideoFolder,sprintf('ICM_%s.avi',datestr(now,'yyyymmdd_HHMMSS'))), 'Uncompressed AVI');
% % diskLogger = VideoWriter(fullfile(handles.cp.VideoFolder,sprintf('ICM_%s.avi',datestr(now,'yyyymmdd_HHMMSS'))), 'Archival');
% diskLogger.FrameRate = handles.cp.FPS;
% handles.video.DiskLogger = diskLogger;
      
%% set measurement parameters by default
pb_DefaultMeasParameters_Callback(hObject, eventdata, handles);
handles = updateParameters(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Real_Time_ICM wait for user response (see UIRESUME)
% uiwait(handles.Real_Time_ICM);% xyz note: by default it is blocked, but some
        ...website suggest turning it on. Yes, should be on.

% --- Outputs from this function are returned to the command line.
