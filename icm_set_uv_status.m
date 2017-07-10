function icm_set_uv_status(onOrOff,uvLevel)

% result folder
global uv
global g_tImaqStart
global g_uvStatus
global hdvp
global DMDimage

if onOrOff == 1
    %% Prepare DMD as a 2nd monitor to display the upcoming image
    % Get monitors (including DMD which is also a monitor) information
    % use get(0,'MonitorPosition') to get the locations and size of the monitors
    % for details http://www.mathworks.com/help/vision/ref/vision.deployablevideoplayer-class.html
    % a = [Left-bottom width, Left-bottom height, screen width, screen height]
    a = get(0,'MonitorPosition');
    hdvp = vision.DeployableVideoPlayer;
    hdvp.Size = 'Full-Screen'; % this command set the video to be displayed in full-screen
    hdvp.Location = a(2,1:2); % this defines the DMD location where the image should display.
    % hdvp.Location = a(1,1:2); % this tests the code on primary monitor
    
    %% Display DMD bitmap 
    step(hdvp,DMDimage);
    
    %% UV lighting and exposure finally begins here
    uv = UVConn('COM3');
    UVSetIrisLevel(uv, uvLevel);
    UVShutterOpen(uv);
    g_tImaqStart = tic;
    g_uvStatus = 1; % flag UV light is on
    disp('open uv')
    
else
    % Close UV Shutter
    UVShutterClose(uv);
    
    % flag UV light is off again after being "on", not "0" so that acquisition and measurement could go on to capture dark curing
    g_uvStatus = 2; 
    disp('close uv');
    
    % Disconnect UV
    UVDisc(uv);

    % update the Interferogram Status Bar for info
%     set(handles.st_InterferogramStatusBar,'String','ON Target. Controller closed UV shutter. Stop measurement when ready.');
    
end    
