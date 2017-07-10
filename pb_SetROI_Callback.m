function pb_SetROI_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SetROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global g_POI

%% Change button value 
%     set(handles.pb_SetROI,'String','Stop & Save Measurement');
%     set(handles.pb_AcquireAVI,'Enable','off');

g_POI = []; % POI (Points of Interest): region to be measured

%% Set ROI and start measuring..            
%% Method 1: NO user interface, Codes specify certain points
%... centerpoint (h,w)=(220,295)estimated on 3/15/2016 setup
% %             hc = 215;
% %             wc = 295;
% single pixel
hc = 210;
wc = 285;
g_POI = [g_POI; [hc wc]];

% % three pixels
% g_POI = [200 275; 210 285; 220 295];

%             % multiple points of horizontal line
%             wc = [280:5:295]';
%             hc = 240*ones(length(wc),1);
%             g_POI = [g_POI; combvec([hc,wc])];

%             % multiple points of vertical line
%             hc = [140:5:280]';
%             wc = 305*ones(length(hc),1);
%             g_POI = [g_POI; combvec([hc,wc])];

% multiple points of Rectangle Area
%             wc = [280:5:295];
%             hc = [225:5:240];
%             g_POI = [g_POI; combvec(hc,wc)'];

%             % to filter neighboring pixels in the ROI for denoising
%             hc = 220;
%             wc = 295;
%             ROI = combvec([hc-1 hc hc+1],[wc-1 wc wc+1]);
%             g_POI = [g_POI; ROI'];
%         
%         %% Method 2: User Interface Select ROI of Rectangle
%             %---- User Input: select a rectangle on figure
%             handles.rect = round(getrect(handles.Interferogram));        
% 
%             % update the selected area to be measured
%             handles.toggle_POI_center = get(handles.checkbox_CenterPointROI,'Value');
%             handles.toggle_POI_corner = get(handles.checkbox_CornerPointsROI, 'Value');
%             handles.toggle_POI_horizon_line = get(handles.checkbox_hCenterlineROI, 'Value');
%             handles.toggle_POI_vertical_line = get(handles.checkbox_vCenterlineROI, 'Value');
% 
%             %ws:width start; we:width end; hs:height start; he:height end
%             ws = handles.rect(1);
%             we = ws + handles.rect(3)-1;
%             hs = handles.rect(2);
%             he = hs + handles.rect(4)-1;
%             % center point
%             wc = round(0.5*(ws + we));
%             hc = round(0.5*(hs + he));
% 
%             guidata(hObject, handles);
% 
%             if handles.toggle_POI_center == 1            
%                 g_POI = [g_POI; [hc wc] ];
%             end
% 
%             if handles.toggle_POI_corner == 1
%                 g_POI = [g_POI; [hs ws]; [hs we]; [he ws]; [he we] ];
%             end
% 
%             if handles.toggle_POI_horizon_line == 1
%                 n = handles.rect(3);
%                 g_POI = [g_POI; [repmat(hc, n, 1) (ws : we)'] ];
%             end
% 
%             if handles.toggle_POI_vertical_line == 1
%                 n = handles.rect(4);
%                 g_POI = [g_POI; [(hs : he)' repmat(wc, n, 1)] ];
%             end
% 
%             %----- display ROI in the figure
%             child=get(gca, 'Children');
%             if size (child,1)==1
%             %     rectangle('Position', handles.shape, 'Curvature', [1 1], 'EdgeColor', 'b');
%                 rectangle('Position', handles.rect,...
%                     'Curvature',[0, 0],...
%                     'EdgeColor', 'r',...
%                     'LineWidth', 3,...
%                     'LineStyle','-')
%             else
%                 set(child(1),'Position', handles.rect);
%             end
% 
%             % hold on;
%             % plot(g_POI', 'r.','MarkerSize',2) 
% 
%             guidata(hObject, handles);

%% returns g_POI with all points of interest to be measured
... 2-by-nPOI matrix: Row 1 is height and Row 2 is width coordinations
g_POI = transpose(g_POI);
disp(g_POI);

if get(handles.ppmenu_RealTimeOrOffline, 'Value') == 1 % Real-time
    set(handles.pb_AcquireAVI,'Enable','on');
    set(handles.st_InterferogramStatusBar,'String','Next Steps: 1.Generate DMD Bitmap;2.Setup control;3.Click "Acquire & Analyze".');
elseif get(handles.ppmenu_RealTimeOrOffline, 'Value') == 2 % Offline
    handles.mask = 1; % flag ROI is set
    global MeasBeginFrame
    MeasBeginFrame = handles.CFrameInd;
    global foiIdx
    foiIdx = 0;
    global RunNo
    RunNo = 0;
    set(handles.st_InterferogramStatusBar,'String','Next Steps: Click "Play" to start measurement of ROI pixels');
end

guidata(hObject, handles);

