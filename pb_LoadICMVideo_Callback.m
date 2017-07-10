%% --- Executes on button press in pb_LoadICMVideo.
function pb_LoadICMVideo_Callback(hObject, eventdata, handles)
% hObject    handle to pb_LoadICMVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clear memory for new data
global ICMvidData
ICMvidData = []; % clear previously loaded ICMvidData in workspace

% get user file
[handles.AVIfilename, handles.AVIpathname] = uigetfile('*.avi', 'Select an ICM Video file to open');
if handles.AVIfilename == 0, return,end; 

global imageTime % Time stamp of acquired video frames
global frameIdx_uvOpen % Frame Index when UV opened
global frameIdx_uvClose % Frame Index when UV closed

[handles.IMAQfilename, handles.IMAQpathname] = uigetfile('*.mat', 'Load the image acquistion file:RT_ECPL_ICM_*.mat file');
if handles.IMAQfilename == 0, return,end; 
load(fullfile(handles.IMAQpathname,handles.IMAQfilename));


% % for debug on Mar-21-2016 only to get "frameIdx_uvOpen" and
% % "frameIdx_uvClose"
% [handles.ECPLfilename, handles.ECPLpathname] = uigetfile('*.mat', 'Load the ECPL file:*ECPL*.mat file');
% if handles.ECPLfilename == 0, return,end; 
% load(fullfile(handles.ECPLpathname,handles.ECPLfilename));


% Read the selected video data
handles.AVIvideoPath = fullfile(handles.AVIpathname, handles.AVIfilename);
[AVIpathstr,AVIname,ext] = fileparts(handles.AVIvideoPath);
handles.AVIdataPath = fullfile(AVIpathstr,[AVIname,'.mat']); % video data

% If: Video was loaded and read before, no need to extract the massive video data again
if exist(handles.AVIdataPath,'file') 
    load(handles.AVIdataPath,'ICMvidData');

% else: First time to read the AVI file and save its data for future uses to avoid repeated consumption of computation resource
else 
    ICMvideoReader = VideoReader(handles.AVIvideoPath);

    %% Method 1: use new function "readFrame"
    FrameIdx = 1;
    while hasFrame(ICMvideoReader)
        vidFrame = readFrame(ICMvideoReader);
        v3dSrc(:,:,FrameIdx) = vidFrame;
        FrameIdx = FrameIdx+1;
    end
    ICMvidData=uint8(v3dSrc);
    
%     %% Alternative method: use old function "read"
%     v3dSrc = ICMvideoReader.read();
%     % Remove gray scale dim
%     if size(v3dSrc,3)==1 %Grayscale AVI: height-by-width-by-1-by-frames
%         ICMvidData = squeeze(v3dSrc);
%     elseif size(v3dSrc,3)==3 %4-D RGB images:height-by-width-by-3-by-frames
%         ICMvidData = squeeze(v3dSrc(:,:,1,:));
%     end

    % save AVI data to file
    save(handles.AVIdataPath, 'ICMvidData');
end

handles.AVIheight = size(ICMvidData,1); % height of video
handles.AVIwidth = size(ICMvidData,2); % width of video
handles.AVInFrame = size(ICMvidData,3); % # of frames in the video

% update the Interferogram Status Bar for info
set(handles.st_InterferogramStatusBar,'String',sprintf('Loaded %s.avi & .mat files at %s ', AVIname,AVIpathstr));
% display the first frame as initial
handles.CFrameInd = 1; % current frame index
CFrame=ICMvidData(:,:,handles.CFrameInd);
axes(handles.Interferogram);
imshow(CFrame);
% set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'xtick', [], 'ytick', []);
set(handles.st_DisplayedFrame,'String',sprintf('Frame %d of %d',handles.CFrameInd,handles.AVInFrame));

handles.mask = 0;
set(handles.st_InterferogramStatusBar,'String','To measure offline, click "Set ROI" then "Play".   To replay only, click "Play"');

guidata(hObject, handles);
pause(0.5);