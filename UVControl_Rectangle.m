function UVControl_Rectangle(recWidth,recHeight,iris,grayscale,ExpTime)
%% This function projects a rectangle with recWidth and recHeight [px] for exposure time [s].
... e.g. UVControl_Rectangle(250,250,22,255,10)
% recWidth,recHeight: width and height of the rectangle (e.g. 250 or 500)
% Iris level between 1 and 95 [%]. 
%GrayScale: grayscale of the bitmap (255=black=DMD on, 0=white=DMD off)
% UV Intensity Control: to test the intesntiy effects, we may set grayscale = 255 (black), and
...vary iris levels.

%% file directory
rootDir = pwd;
DMD_Folder = fullfile(rootDir, 'DMD'); 
% create a folder if it doesn't exist
if ~exist(DMD_Folder, 'dir')
    mkdir(rootDir, 'DMD');
end

%% Get monitors (including DMD which is also a monitor) information
 % a = [Left-bottom width, Left-bottom height, screen width, screen height]
a = get(0,'MonitorPosition');
hdvp = vision.DeployableVideoPlayer;
hdvp.Size = 'Full-Screen'; % this command set the video to be displayed in full-screen
% hdvp.Location = a(2,1:2); % this defines the DMD location where the image should display.
hdvp.Location = a(1,1:2); % this tests the code on primary monitor
% use get(0,'MonitorPosition') to get the locations and size of the monitors
% e.g.: get(0,'MonitorPosition')
% 
% ans =
% 
%        -1279         150        1280        1024
%            1           1        1680        1050
% then setting hdvp.Location = [1,1] will get the viodeo to show on the
% second screen; setting hdvp.Location = [-1279,150] will get the video to
% show on the first screen
% see http://www.mathworks.com/help/vision/ref/vision.deployablevideoplayer-class.html
% for details
%% generate a rectangle
    w = 1024; %DMD width [px]
    h = 768; %DMD height [px]
% invert:normally Black is DMD "on", white is DMD "off",
... but in the grey bitmap,greyscale in bitmap is 0 = black, 255 = white
... so, we need to invert the greyscale so that the greater intensity should be greater greyscale
    invert = 1; 

    img = generateRectangle(w,h,recWidth,recHeight,grayscale,invert);
%     generateRectangleVideo(1);

%% Control UV lighting time
    t = 0;
%%--- notes: when connected to Lab computer, 
...uncomment all the following lines with UV control
%     uv = UVConn('COM3');
%     UVSetIrisLevel(uv,iris);
    tic
%     UVShutterOpen(uv);
    while t < ExpTime;
        step(hdvp,img);
        t = toc;
    end
%     UVShutterClose(uv);
%     UVDisc(uv);

end