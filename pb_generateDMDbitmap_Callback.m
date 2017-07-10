% --- Executes on button press in pb_generateDMDbitmap.
function pb_generateDMDbitmap_Callback(hObject, eventdata, handles)
%% This function:
... 1) generates a bitmap by funciton"generateRectangle.m"
    % recWidth,recHeight: width and height of the rectangle (unit: pixels)
    % GrayScale: grayscale of the bitmap (255=black=DMD on, 0=white=DMD off)
... 2) Saves the DMD image to result folder
    
% result folder
global DMDimage

%% Generate a bitmap
    w = 1024; %DMD width [px]
    h = 768; %DMD height [px]
    % invert:normally Black is DMD "on", white is DMD "off",
    ... but in the grey bitmap,greyscale in bitmap is 0 = black, 255 = white
    ... so, we need to invert the greyscale so that the greater intensity should be greater greyscale
    invert = 1; 
    
    handles.recWidth = str2double(get(handles.ed_recWidth,'String'));
    handles.recHeight = str2double(get(handles.ed_recHeight,'String'));
    handles.recGrayscale = str2double(get(handles.ed_recGrayscale,'String'));

    % generate and save a DMD bitmap
    DMDimage = generateRectangle(w,h,handles.recWidth,handles.recHeight,handles.recGrayscale,invert);
%         DMDimage = generateRectangle(1024,768,350,250,255,1);

%% save image file
% % file directory
% rootDir = pwd;
% DMD_Folder = fullfile(rootDir, 'DMD'); 
DMD_Folder = fullfile(handles.cp.ResultFolder, 'DMD'); 
% create a folder if it doesn't exist
if ~exist(DMD_Folder, 'dir')
    mkdir(DMD_Folder);
end

imwrite(DMDimage,fullfile(DMD_Folder,strcat('DMD_Bitmap_Rec_',datestr(now,'yyyymmdd_HHMMSS'),'.png')));
recWidth = handles.recWidth;
recHeight = handles.recHeight;
recGrayscale = handles.recGrayscale;
save(strcat(handles.cp.ResultFolder,strcat('\DMD_Bitmap_Rec_',datestr(now,'yyyymmdd_HHMMSS'),'.mat')),'DMDimage',...
         'recWidth','recHeight','recGrayscale'); 

%% update guidata
guidata(hObject, handles);
end

