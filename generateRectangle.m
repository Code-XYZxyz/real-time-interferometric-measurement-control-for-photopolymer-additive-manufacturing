function [ img ] = generateRectangle(w,h,recWidth,recHeight,GreyScale,invert)
%% this function generates a grey bitmap with constant greyscales 
% for all pixels
... (h,w):bitmap height and width, 
    ...w = 1024; %DMD width
    ...h = 768; %DMD height
... recWidth,recHeight: width and height of the rectangle (e.g. 250 or 500)
... GrayScale: grayscale of the bitmap
... invert:normally Black is DMD "on", white is DMD "off",
... but in the grey bitmap,greyscale in bitmap is 0 = black, 255 = white
... so, we need to invert the greyscale so that the greater intensity should be greater greyscale

%%
if nargin < 6
    invert = 1; % by default to invert it to match the greyscale with the intensity scale
end
    img = uint8(zeros(h,w,1));
    img(:) = 255;
    hBorder = round((h-recHeight)/2);
    wBorder = round((w-recWidth)/2);
%% Greyscale Bitmap generation
if invert
    % 255 is black, DMD max intensity; 0 is white, DMD main intensity
    img(hBorder:hBorder+recHeight-1,wBorder:wBorder+recWidth-1) = 255-GreyScale; 
else
    img(hBorder:h-hBorder,wBorder:w-wBorder) = GreyScale;
end
% %% save image file
% % % file directory
% % rootDir = pwd;
% % DMD_Folder = fullfile(rootDir, 'DMD'); 
% DMD_Folder = fullfile(ResultFolder, 'DMD'); 
% % create a folder if it doesn't exist
% if ~exist(DMD_Folder, 'dir')
%     mkdir(DMD_Folder);
% end
% 
% imwrite(img,fullfile(DMD_Folder,strcat('DMD_Bitmap_Rec_',datestr(now,'yyyymmdd_HHMMSS'),'.png')));
end