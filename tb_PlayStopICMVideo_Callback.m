%% --- Executes on button press in tb_PlayStopICMVideo.
function tb_PlayStopICMVideo_Callback(hObject, eventdata, handles)
% hObject    handle to tb_PlayStopICMVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hImage
global ICMvidData % all data in the loaded ICM video
global MeasBeginFrame % index of first frame in the measurement span of video data
global foiIdx % foi (frame of interest) Index
% global dataY % time sequence of grayscales
global RunNo % number of runs of rolling fit

% Loaded data at "pb_LoadICMVideo_Callback.m"
global imageTime % Image time stamp file loaded at "pb_LoadICMVideo_Callback.m"
global frameIdx_uvOpen % Frame Index when UV opened
global frameIdx_uvClose % Frame Index when UV closed

global foiTimeAbs_dummy % foi time with missing frames times make-up

global g_POI % Pixels of interest returned by "pb_SetROI_Callback.m"

%==== Measurement result return: rolling fit coeff, online height estimate
% nPOI-by-1 structure array of all points measurement
% Each point structure has fields: 'PixelHeightWidth', 'rawY','fitY','FittedCoeffs'
... 'CureFlags', 'Idx_FailFit', 'Times', 'Heights', 'Freq_w','Freq','FreqCumSum'
% Initialized at callback function "pb_SetROI.m"
global MeasureRet
global RunNo_uvClose 
global zExposedNorminal % ICM Measured Average Cured Heights across ROI when UV closes
global zDarkNorminal % ICM Measured Average Cured Heights across ROI after UV closes at the end of acquisition
global NumThresholdPixels_Array
global NumDarkPixels_Array

if isempty(ICMvidData),return, end; % if no offline ICM video loaded

handles = updateParameters(handles);
MeasStruct = icm_init_measure_ret(handles.cp);
nPOI = size(g_POI, 2);
MeasureRet = repmat(MeasStruct,nPOI,1); 

% get(hObject,'Value')= 1, on (video playing); = 0, off (video stops)
%% if "play" is pressed, toggle button turns on
while get(hObject,'Value') && (handles.CFrameInd < handles.AVInFrame + 1)
    set(handles.tb_PlayStopICMVideo, 'String', 'Stop','Enable', 'on');
    set(handles.st_InterferogramStatusBar,'String','Replaying ICM&M video');
    set(handles.pb_ICMVideo1stFrame, 'Enable', 'off');
    set(handles.pb_ICMVideoPrevFrame, 'Enable', 'off');
    set(handles.pb_ICMVideoNextFrame, 'Enable', 'off');
    set(handles.pb_ICMVideoEndFrame, 'Enable', 'off');
    handles = updateParameters(handles);
    %% play acquired ICM video
    CFrame=ICMvidData(:,:,handles.CFrameInd);
    if handles.CFrameInd == 1
        hImage = imshow(CFrame, 'Parent', handles.Interferogram);
    else
        set(hImage, 'CData', CFrame);
    end
    set(handles.st_DisplayedFrame,'String',sprintf('Frame %d of %d',handles.CFrameInd,handles.AVInFrame));
    handles.CFrameInd = handles.CFrameInd+1; % prepare to show next frame,need reduce by "1" when stop to play
    guidata(hObject, handles);
    
    %% offline analysis of the video data starts when ROI is set
   if handles.mask == 1
        %     if strcmp(get(handles.pb_SetROI,'String'),'Stop & Save Measurement')
        set(handles.tb_PlayStopICMVideo, 'Enable', 'off'); % should stop measurement before stop video
        set(handles.st_InterferogramStatusBar,'String','Offline ICM&M measurement and analysis is ongoing');
        foiIdx = foiIdx + 1;
        
        %% nPOI: number of Points of Interest in a frame
        ...global g_POI is returned by "pb_SetROI_Callback.m"
        ...global g_POI: all points of interest in a frame ROI
           nPOI = size(g_POI, 2);
       
%        %% --- make-up dataX (image time) with missing frames, will be used for predicting dataY (grayscales)
%                  % if elasped time since last frame longer than 2 IMAQ periods, need make up for missing data
%                     if (RunNo > 2) && (imageTime(foiIdx)-imageTime(foiIdx-1)> 5/handles.cp.FPS)
%                         foiTimeAbs_miss = imageTime(foiIdx-1)+ 1/handles.cp.FPS: (1/handles.cp.FPS):imageTime(foiIdx)-1/handles.cp.FPS;
%                         foiTimeAbs_dummy = [foiTimeAbs_dummy; foiTimeAbs_miss';imageTime(foiIdx)]; 
%                         dataX_miss = foiTimeAbs_miss - imageTime(MeasBeginFrame);
%                      else
%                         foiTimeAbs_dummy = [foiTimeAbs_dummy;imageTime(foiIdx)];                      
%                     end
%           %%----------
       
        %% extract grayscale data for the POI using median filter to denoise
        for iPoint = 1:nPOI
            % single point identification by coordinations (height, width)
            h = g_POI(1, iPoint); % height coordinate
            w = g_POI(2, iPoint); % width cooridnate
            
            % time series of intensity data
            % Note: ICMvidData is a matrix of Height-by-Width-by-FrameNum
%             area = ICMvidData((h-3):(h+3), (w-3):(w+3), foiIdx+MeasBeginFrame-1); % 7X7 filter
            area = ICMvidData((h-2):(h+2), (w-2):(w+2), foiIdx+MeasBeginFrame-1); % 5X5 filter
%             area = ICMvidData(h, w, foiIdx+MeasBeginFrame-1); % 1X1 filter
            dataY_foi = double(median(area(:)));
            MeasureRet(iPoint).PixelHeightWidth = [h;w];
            MeasureRet(iPoint).rawY = [MeasureRet(iPoint).rawY; double(ICMvidData(h, w,foiIdx+MeasBeginFrame-1))]; % pixel grayscale
            MeasureRet(iPoint).dataY = [MeasureRet(iPoint).dataY; dataY_foi]; % filtered grayscale

%             % time series of intensity data
% %             %OLD method: single pixel intensity
% %             dataY = double(squeeze(ICMvidData(h, w, foiIdx)));
% %             MeasureRet(iPoint).rawY = dataY;
% 
%             %NEW method: median filtering neighboring pixels for single
%             ... pixel intensity to reduce noise
% %             % 11X11 filter
% %             dataY_foi = double(median(median(ICMvidData(h-5:h+5, w-5:w+5,foiIdx+MeasBeginFrame-1))));
% 
% %             % 9X9 filter
% %             dataY_foi = double(median(median(ICMvidData(h-4:h+4, w-4:w+4,foiIdx+MeasBeginFrame-1))));
% 
%             % 7X7 filter
%             dataY_foi = double(median(median(ICMvidData(h-3:h+3, w-3:w+3,foiIdx+MeasBeginFrame-1))));
% %             
% %             % 5X5 filter
% %             dataY_foi = double(median(median(ICMvidData(h-2:h+2, w-2:w+2,foiIdx+MeasBeginFrame-1)))); 
% 
% %             % 3X3 filter
% %             dataY_foi = double(median(median(ICMvidData(h-1:h+1, w-1:w+1,foiIdx+MeasBeginFrame-1)))); 
% 
% %             % 1X1 filter
% %             dataY_foi = double(median(median(ICMvidData(h-0:h+0, w-0:w+0,foiIdx+MeasBeginFrame-1)))); 
% 
%             MeasureRet(iPoint).rawY = [MeasureRet(iPoint).rawY;dataY_foi];
            
%              %%--- make-up dataY with missing frames
%                  % if elasped time since last frame longer than 2 IMAQ periods, need make up for missing data
%                     if (RunNo > 2) && (imageTime(foiIdx)-imageTime(foiIdx-1)> 5/handles.cp.FPS)
%                         if isfield(MeasureRet(iPoint).lastFitRet,'fitobject')
%                             if (MeasureRet(iPoint).lastFitRet.fitStatus ~=0 && MeasureRet(iPoint).lastFitRet.fitStatus ~=3 && MeasureRet(iPoint).lastFitRet.fitStatus ~=40)...
%                                     ||(MeasureRet(iPoint).lastFitRet.fitStatus ==0 && MeasureRet(iPoint).lastFitRet.freq~=0 ) % small freq (<0.1Hz)in threshold acceptable
%                                 PredictY = feval(MeasureRet(iPoint).lastFitRet.fitobject, dataX_miss);
%                                 MeasureRet(iPoint).dummyY = [MeasureRet(iPoint).dummyY; PredictY;dataY_foi];
%                             else
%                                PredictY = MeasureRet(iPoint).lastFitRet.I0*ones(length(dataX_miss),1);
%                                 MeasureRet(iPoint).dummyY = [MeasureRet(iPoint).dummyY; PredictY;dataY_foi]; 
%                             end
%                         else
%                                 PredictY = MeasureRet(iPoint).lastFitRet.I0*ones(length(dataX_miss),1);
%                                 MeasureRet(iPoint).dummyY = [MeasureRet(iPoint).dummyY; PredictY;dataY_foi]; 
%                         end
%                     else
%                         MeasureRet(iPoint).dummyY = [MeasureRet(iPoint).dummyY;dataY_foi];
%                     end
%              %%----------
            
        end
         %% Rolling fit: when sufficient samples
            ... and when every measurement period arrived
                % handles.cp.MeasPeriodSamples: online update the model parameters every "MeasPeriodSamples" new data is
            ... acquired, and meanwhile predict next set of "MeasPeriodSamples" data.
                % handles.cp.MeasPeriodSamples = str2double(get(handles.ed_MeasPeriodSamples,'String'));
            %             if (foiIdx >= handles.cp.MovingHorizonL)&&(mod(foiIdx-handles.cp.MovingHorizonL, handles.cp.MeasPeriodSamples) == 0)
            if ((foiIdx >= handles.cp.SamplesNumB4Measure) && (mod(foiIdx-handles.cp.SamplesNumB4Measure, handles.cp.MeasPeriodSamples) == 0))...
                    || (foiIdx == handles.AVInFrame)


%                 % nPOI: number of Points of Interest in a frame
%                 ...global g_POI is returned by "pb_SetROI_Callback.m"
%                 ...global g_POI: all points of interest in a frame ROI
%                 nPOI = size(global g_POI, 2);

                %% point-by-poiont analysis
                RunNo = RunNo + 1; % Run number of rolling fit&prediction
                dataX = imageTime(MeasBeginFrame : foiIdx+MeasBeginFrame-1)-imageTime(MeasBeginFrame);
%                 dataX = imageTime(MeasBeginFrame : MeasBeginFrame + foiIdx-1)-imageTime(MeasBeginFrame);
%                 dataX_dummy = foiTimeAbs_dummy(1:end)-foiTimeAbs_dummy(1); % time of each foi relative to start of measurement with missing data make-up
%                 dataX = imageTime(end-foiIdx+1:end)-imageTime(end-foiIdx); % time of each foi relative to start of measurement
%                 dataX = foiTimeAbs(2:end)-foiTimeAbs(1); % time of each foi relative to start of measurement
                % dataX = (1:foiIdx)' / handles.cp.FPS; % time "t" assume constant FPS
                % Curve fitting parameters
                params.rSquare = handles.cp.GOF_rSquare;
                params.FPS = handles.cp.FPS;
                params.halfLife = handles.cp.HalfLife;
                params.MHL = handles.cp.MovingHorizonL;
                params.MeasPeriod = handles.cp.MeasPeriodSamples;
                params.uvIris = handles.cp.uvIris;
                params.f_max = handles.cp.f_max;
                params.f_diff_max = handles.cp.f_diff_max;
                params.RunNo= RunNo;
                params.frameIdx   = foiIdx+MeasBeginFrame-1;

                % Mark the end of exposed curing, i.e., start of dark curing 
                if  (foiIdx+MeasBeginFrame-1 >= frameIdx_uvClose) && (isempty(RunNo_uvClose))
                    RunNo_uvClose = RunNo;% flag the frame number in acquired video when UV closes
                end
                params.RunNo_uvClose = RunNo_uvClose;
                 
%                %%============ Begin: Serial Computation ============
%                 % Start curve fitting and heights calculation
%                 for iPoint = 1:nPOI
%                     % single point identification by coordinations (height, width)
%                     h = g_POI(1, iPoint); % height coordinate
%                     w = g_POI(2, iPoint); % width cooridnate
%                     MeasureRet(iPoint).PixelHeightWidth = [h;w];
% 
%                     % time series of intensity data
% %                     dataY = double(squeeze(ICMvidData(h, w,  MeasBeginFrame:foiIdx+MeasBeginFrame-1)));
%                     dataY = MeasureRet(iPoint).rawY; % without missing data imputation
%                     dataY_dummy = MeasureRet(iPoint).dummyY; % WITH missing data imputation
%                     MeasureRet(iPoint).dataX = dataX;
%                     MeasureRet(iPoint).dataX_dummy = dataX_dummy;
% 
%                   %% Curve Fitting
%                     % Rolling fit with "fourier1" returns 4 coefficients y=a0+a1*cos(px)+b1*sin(px)
%                     % fitRollRet = icmFit(trainX, trainY, trainW);
%                     if RunNo == 1
%                         prevFitRet.fitStatus = 0;
%                         prevFitRet.time = 0;
%                         prevFitRet.freq = 0;
%                         prevFitRet.frameIdx = 0;
%                         prevFitRet.foiIdx_dummy = 0;
%                         prevFitRet.firstValidFoiIdx = 0;
%                         prevFitRet.freq = 0;
%                         prevFitRet.movingHorizon = handles.cp.MovingHorizonL;
%                         prevFitRet.halfLife = handles.cp.HalfLife;
%                     else
%                         prevFitRet = MeasureRet(iPoint).lastFitRet;
%                     end
% %                     fitRollRet = icmFit2(dataX, dataY, params, prevFitRet); % without missing data imputation
%                     fitRollRet = icmFit2(dataX_dummy, dataY_dummy, params, prevFitRet);% WITH missing data imputation
%                     MeasureRet(iPoint).lastFitRet = fitRollRet;
% 
% 
%                     % save the fitting coefficients, i.e., online estimates of parameters
%                     coeffs = [fitRollRet.fitStatus,fitRollRet.fitgof.rsquare,...
%                         fitRollRet.I0, fitRollRet.I1, fitRollRet.freqW, fitRollRet.freq,...
%                         fitRollRet.movingHorizon, fitRollRet.halfLife];
%                     MeasureRet(iPoint).FittedCoeffs = [MeasureRet(iPoint).FittedCoeffs; coeffs];
%                     MeasureRet(iPoint).Freq_w = MeasureRet(iPoint).FittedCoeffs(:,5);
%                     MeasureRet(iPoint).Freq = MeasureRet(iPoint).FittedCoeffs(:,6);
%                     
% %                     if fitRollRet.fitStatus ~=0 && fitRollRet.fitStatus ~=3
%                     if isfield(fitRollRet,'fitobject')
%                         if (fitRollRet.fitStatus ~=0 && fitRollRet.fitStatus ~=3 && fitRollRet.fitStatus ~=40)...
%                                 ||(fitRollRet.fitStatus ==0 && fitRollRet.freq~=0 ) % small freq (<0.1Hz)in threshold acceptable
%                             newFitY = feval(fitRollRet.fitobject, dataX( (prevFitRet.frameIdx+1):fitRollRet.frameIdx));
% %                             newFitY = feval(fitRollRet.fitobject, dataX( end-handles.cp.MeasPeriodSamples+1:end));
%                             MeasureRet(iPoint).fitY = [MeasureRet(iPoint).fitY; newFitY];
%                         else
%                            newFitY = fitRollRet.I0*ones(fitRollRet.frameIdx-prevFitRet.frameIdx,1);
%                             MeasureRet(iPoint).fitY = [MeasureRet(iPoint).fitY; newFitY]; 
%                         end
%                     else
%                         newFitY = fitRollRet.I0*ones(fitRollRet.frameIdx-prevFitRet.frameIdx,1);
%                         MeasureRet(iPoint).fitY = [MeasureRet(iPoint).fitY; newFitY];
%                     end
% 
%                      % Mark the start of curing, i.e., the end of threshold
%                     if  (MeasureRet(iPoint).CureFlags.CureFlag_RunNo==0) && (fitRollRet.firstValidFoiIdx~=0)
%                         MeasureRet(iPoint).CureFlags.CureFlag_RunNo = RunNo;
%                         MeasureRet(iPoint).CureFlags.CureFlag_FrameIdx = fitRollRet.firstValidFoiIdx;
%                     end
%                    % Flag the runs of failed curve fitting, which has low R-square
%                     ...and may yield frequency outlier
%                     if fitRollRet.fitStatus ~= 2
%                         MeasureRet(iPoint).Idx_FailFit = [MeasureRet(iPoint).Idx_FailFit; RunNo];
%                     end
% 
%                  %% Height Estimation: growth by integration
%                    T_Int = dataX(end)- prevFitRet.time;
%                    % array of measurement time(s) per point, RunNo-by-1 matrix
%                    MeasureRet(iPoint).Times = [MeasureRet(iPoint).Times;dataX(end)];
%                    
%                    % phase(unit: 2Pi: time cumulative sum of frequency-by-time for height estimation
% %                    %--- method 1: mixed use of trapzoidal and local value
% %                    if T_Int < 1.5*handles.cp.MeasPeriodSamples/handles.cp.FPS
% %                        MeasureRet(iPoint).Phase2Pi = MeasureRet(iPoint).Phase2Pi+T_Int*fitRollRet.freq;
% %                    else  % if too long interval, use average freq 
% %                        MeasureRet(iPoint).Phase2Pi = MeasureRet(iPoint).Phase2Pi+T_Int*(fitRollRet.freq+prevFitRet.freq)/2;
% %                    end
% %                    %--- method 2: always use trapzoidal rule, i.e., midpoint
% %                    MeasureRet(iPoint).Phase2Pi = MeasureRet(iPoint).Phase2Pi+T_Int*(fitRollRet.freq+prevFitRet.freq)/2;
%                    %--- method 3: always use local frequency only
%                    MeasureRet(iPoint).Phase2Pi = MeasureRet(iPoint).Phase2Pi+T_Int*fitRollRet.freq;
%                    
%                    % array of cured heights
%                    z = handles.cp.Wavelength/(2*(handles.cp.n_m-handles.cp.n_L))*MeasureRet(iPoint).Phase2Pi;
%                    MeasureRet(iPoint).Heights = [MeasureRet(iPoint).Heights;z];
%                  
% %                  %--- old method: assume constant FPS
% %                     T_Int = handles.cp.MeasPeriodSamples/handles.cp.FPS;
% %                     % array of run time per point, RunNo-by-1 matrix
% %                     MeasureRet(iPoint).Times = [MeasureRet(iPoint).Times;(RunNo-1)*T_Int];
% %                     % cumulative sum of frequencies for height estimation
% %                     MeasureRet(iPoint).FreqCumSum = sum(MeasureRet(iPoint).Freq);
% %                     % array of cured heights
% %                     z = handles.cp.Wavelength*T_Int/(2*(handles.cp.n_m-handles.cp.n_L))*MeasureRet(iPoint).FreqCumSum;
% %                     MeasureRet(iPoint).Heights = [MeasureRet(iPoint).Heights;z];
%                  
% 
%                 end
%              %%============ END: Serial Computation ============

% % -------------------------------------------Dividing Line----------------------------------------------------------------------

            %% parallel computing for multi-pixel measurement
            POI = g_POI;
            MHL = handles.cp.MovingHorizonL;
            HalfLife = handles.cp.HalfLife;
            MeasPeriodSamples=handles.cp.MeasPeriodSamples;
            Wavelength = handles.cp.Wavelength;
            n_m = handles.cp.n_m;
            n_L = handles.cp.n_L;
        
%         % before04-08-2016
%         for iPoint = 1:nPOI
%             prevFitRet(iPoint) = MeasureRet(iPoint).lastFitRet;
%         end

        % Updated on 04-08-2016: majority voting to decide curing start
        NumThresholdPixels = 0; % Number of Pixels that are in threshold
        NumDarkPixels = 0; % Number of Pixels that are in threshold
        for iPoint = 1:nPOI
            prevFitRet(iPoint) = MeasureRet(iPoint).lastFitRet;
            if prevFitRet(iPoint).fitStatus == 0
                NumThresholdPixels = NumThresholdPixels+1;
            elseif prevFitRet(iPoint).fitStatus == 3
                NumDarkPixels = NumDarkPixels+1;
            end
        end
        NumThresholdPixels_Array = [NumThresholdPixels_Array; NumThresholdPixels];
        NumDarkPixels_Array = [NumDarkPixels_Array;NumDarkPixels];
        
        if NumThresholdPixels >= floor(2*nPOI/3) % majority voting for curing start (threshold end)
            for iPoint = 1:nPOI
                if prevFitRet(iPoint).fitStatus ~= 0
                    MeasureRet(iPoint).lastFitRet.fitobject = [];
                    MeasureRet(iPoint).lastFitRet.fitgof.rsquare = 1;
                    MeasureRet(iPoint).lastFitRet.movingHorizon = params.MHL;
                    MeasureRet(iPoint).lastFitRet.halfLife = params.halfLife;                
                    MeasureRet(iPoint).lastFitRet.firstValidFoiIdx = 0; % NOT curing frame yet
                    MeasureRet(iPoint).lastFitRet.fitStatus = 0; % NO fitting yet
                    % retrun DC (Direct Current) values
                    MeasureRet(iPoint).lastFitRet.I0 = mean(MeasureRet(iPoint).dataY((end-MeasPeriodSamples+1):end));
                    MeasureRet(iPoint).lastFitRet.I1 = 0;
                    MeasureRet(iPoint).lastFitRet.freqW = 0;
                    MeasureRet(iPoint).lastFitRet.freq = 0;
                    prevFitRet(iPoint) = MeasureRet(iPoint).lastFitRet;
%                     if RunNo ~= 1
                        MeasureRet(iPoint).Freq_w = zeros(length(MeasureRet(iPoint).Freq_w),1);
                        MeasureRet(iPoint).Freq = zeros(length(MeasureRet(iPoint).Freq),1);
                        MeasureRet(iPoint).Phase2Pi = 0;
                        MeasureRet(iPoint).Heights = zeros(length(MeasureRet(iPoint).Heights),1);
                        MeasureRet(iPoint).zExposed = 0;
                        MeasureRet(iPoint).fitY((end-MeasPeriodSamples+1):end) = MeasureRet(iPoint).lastFitRet.I0*ones(MeasPeriodSamples,1);
%                     end
                end
            end
        end
        if NumDarkPixels >= floor(4*nPOI/5) % majority voting for dark period
            for iPoint = 1:nPOI
                if prevFitRet(iPoint).fitStatus ~= 3 % force the pixels to enter dark curing at the upcoming run
%                     MeasureRet(iPoint).lastFitRet.fitobject = [];
%                     MeasureRet(iPoint).lastFitRet.fitgof.rsquare = 1;
%                     MeasureRet(iPoint).lastFitRet.movingHorizon = params.MHL;
%                     MeasureRet(iPoint).lastFitRet.halfLife = params.halfLife;                
%                     MeasureRet(iPoint).lastFitRet.firstValidFoiIdx = 0; % NOT curing frame yet
                    MeasureRet(iPoint).lastFitRet.fitStatus = 3; % dark curing
%                     % retrun DC (Direct Current) values
%                     MeasureRet(iPoint).lastFitRet.I0 = mean(MeasureRet(iPoint).dataY((end-MeasPeriodSamples+1):end));
%                     MeasureRet(iPoint).lastFitRet.I1 = 0;
%                     MeasureRet(iPoint).lastFitRet.freqW = 0;
%                     MeasureRet(iPoint).lastFitRet.freq = 0;
                    prevFitRet(iPoint) = MeasureRet(iPoint).lastFitRet;
%                     MeasureRet(iPoint).fitY = MeasureRet(iPoint).lastFitRet.I0*ones(MeasPeriodSamples,1);
%                     if ~isempty(MeasureRet(iPoint).Freq_w)
%                         MeasureRet(iPoint).Freq_w(end) = 0;
%                         MeasureRet(iPoint).Freq(end) = 0;
%                         MeasureRet(iPoint).zDark = MeasureRet(iPoint).zDark - (MeasureRet(iPoint).Heights(end) - MeasureRet(iPoint).Heights(end-1));
%                         MeasureRet(iPoint).Heights(end) = MeasureRet(iPoint).Heights(end-1);
%                     end
                end
            end
        end
        
            MeasureRet_par = MeasureRet;
            RunNo_par = RunNo;
        
        % Initialize the local variables for parallel loop
            icmRetStruct=struct('fitobject',[],'fitgof',struct(),'movingHorizon',[],'halfLife',[],...
            'frameIdx_dummy',[],'frameIdx',[],'firstValidFoiIdx',[],'fitStatus',[],'fitoptions',[],...
            'I0',[],'I1',[],'freqW',[],'freq',[],'time',[]);
            fitRollRet = repmat(icmRetStruct,nPOI,1);
            coeffs = cell(nPOI,1);
            newFitY = cell(nPOI,1);
            T_Int = cell(nPOI,1);
            z = cell(nPOI,1);
            
            
            
             % Start curve fitting and heights calculation
%             parpool(4);
                parfor iPoint = 1:nPOI
                    MeasureRet_par(iPoint).dataX = dataX;
%                     MeasureRet_par(iPoint).dataX_dummy = dataX_dummy;

                  %% Curve Fitting
                    % Rolling fit with "fourier1" returns 4 coefficients y=a0+a1*cos(px)+b1*sin(px)
                    fitRollRet(iPoint) = icmFit2(dataX, MeasureRet_par(iPoint).dataY, params, prevFitRet(iPoint)); % with filtered data
             
%                     fitRollRet(iPoint) = icmFit2(dataX, MeasureRet_par(iPoint).rawY, params, prevFitRet(iPoint)); % without missing data imputation
%                     fitRollRet(iPoint) = icmFit2(dataX_dummy, MeasureRet_par(iPoint).dummyY, params, prevFitRet(iPoint));% WITH missing data imputation
                    MeasureRet_par(iPoint).lastFitRet = fitRollRet(iPoint);


                    % save the fitting coefficients, i.e., online estimates of parameters
                    coeffs{iPoint} = [fitRollRet(iPoint).fitStatus,fitRollRet(iPoint).fitgof.rsquare,...
                        fitRollRet(iPoint).I0, fitRollRet(iPoint).I1, fitRollRet(iPoint).freqW, fitRollRet(iPoint).freq,...
                        fitRollRet(iPoint).movingHorizon, fitRollRet(iPoint).halfLife];
                    MeasureRet_par(iPoint).FittedCoeffs = [MeasureRet_par(iPoint).FittedCoeffs; coeffs{iPoint}];
                    MeasureRet_par(iPoint).Freq_w = [MeasureRet_par(iPoint).Freq_w; MeasureRet_par(iPoint).FittedCoeffs(end,5)];
                    MeasureRet_par(iPoint).Freq = [MeasureRet_par(iPoint).Freq; MeasureRet_par(iPoint).FittedCoeffs(end,6)];
                    

                    if ~isempty(fitRollRet(iPoint).fitobject)
                        if (fitRollRet(iPoint).fitStatus ~=0 && fitRollRet(iPoint).fitStatus ~=3 && fitRollRet(iPoint).fitStatus ~=40)...
                                ||(fitRollRet(iPoint).fitStatus ==0 && fitRollRet(iPoint).freq~=0 ) % small freq (<0.1Hz)in threshold acceptable
                            newFitY{iPoint} = feval(fitRollRet(iPoint).fitobject, dataX( (prevFitRet(iPoint).frameIdx+1):fitRollRet(iPoint).frameIdx));
%                             newFitY{iPoint} = feval(fitRollRet(iPoint).fitobject, dataX( end-MeasPeriodSamples+1:end));
                            MeasureRet_par(iPoint).fitY = [MeasureRet_par(iPoint).fitY; newFitY{iPoint}];
                        else
                           newFitY{iPoint} = fitRollRet(iPoint).I0*ones(fitRollRet(iPoint).frameIdx-prevFitRet(iPoint).frameIdx,1);
                            MeasureRet_par(iPoint).fitY = [MeasureRet_par(iPoint).fitY; newFitY{iPoint}]; 
                        end
                    else
                        newFitY{iPoint} = fitRollRet(iPoint).I0*ones(fitRollRet(iPoint).frameIdx-prevFitRet(iPoint).frameIdx,1);
                        MeasureRet_par(iPoint).fitY = [MeasureRet_par(iPoint).fitY; newFitY{iPoint}];
                    end

                     % Mark the start of curing, i.e., the end of threshold
                    if  (MeasureRet_par(iPoint).CureFlags.CureFlag_RunNo==0) && (fitRollRet(iPoint).firstValidFoiIdx~=0)
                        MeasureRet_par(iPoint).CureFlags.CureFlag_RunNo = RunNo_par;
                        MeasureRet_par(iPoint).CureFlags.CureFlag_FrameIdx = fitRollRet(iPoint).firstValidFoiIdx;
                    end
                   % Flag the runs of failed curve fitting, which has low R-square
                    ...and may yield frequency outlier
                    if fitRollRet(iPoint).fitStatus ~= 2
                        MeasureRet_par(iPoint).Idx_FailFit = [MeasureRet_par(iPoint).Idx_FailFit; RunNo_par];
                    end

                 %% Height Estimation: growth by integration
                   T_Int{iPoint} = dataX(end)- prevFitRet(iPoint).time;
                   % array of measurement time(s) per point, RunNo-by-1 matrix
                   MeasureRet_par(iPoint).Times = [MeasureRet_par(iPoint).Times;dataX(end)];
                   
                   % phase(unit: 2Pi: time cumulative sum of frequency-by-time for height estimation
%                    %--- method 1: mixed use of trapzoidal and local value
%                    if T_Int < 1.5*handles.cp.MeasPeriodSamples/handles.cp.FPS
%                        MeasureRet_par(iPoint).Phase2Pi = MeasureRet_par(iPoint).Phase2Pi+T_Int{iPoint}*fitRollRet.freq;
%                    else  % if too long interval, use average freq 
%                        MeasureRet_par(iPoint).Phase2Pi = MeasureRet_par(iPoint).Phase2Pi+T_Int{iPoint}*(fitRollRet.freq+prevFitRet.freq)/2;
%                    end
                   %--- method 2: always use trapzoidal rule, i.e., midpoint
                   MeasureRet_par(iPoint).Phase2Pi = MeasureRet_par(iPoint).Phase2Pi+T_Int{iPoint}*(fitRollRet(iPoint).freq+prevFitRet(iPoint).freq)/2;
%                    %--- method 3: always use local frequency only
%                    MeasureRet_par(iPoint).Phase2Pi = MeasureRet_par(iPoint).Phase2Pi+T_Int{iPoint}*fitRollRet(iPoint).freq;
                   
                   
                  % array of cured heights
%                    % Before 08/24/2016: use constant refractive index
%                    z{iPoint} = Wavelength/(2*(n_m-n_L))*MeasureRet_par(iPoint).Phase2Pi; % constant refractive index

                    % Created on 08/24/2016: use evolving refractive index
                   n_m_evolve = 0.00041*(MeasureRet_par(iPoint).Phase2Pi)+1.49191;% 08/24/2016:calculate evolving refractive index with the model in thesis
                   z{iPoint} = Wavelength/(2*(n_m_evolve-n_L))*MeasureRet_par(iPoint).Phase2Pi; % 08/24/2016: evolving refractive index
                   % end of 08/24 updates
                   
                   MeasureRet_par(iPoint).Heights = [MeasureRet_par(iPoint).Heights;z{iPoint}];
                 
%                  %--- old method: assume constant FPS
%                     T_Int = handles.cp.MeasPeriodSamples/handles.cp.FPS;
%                     % array of run time per point, RunNo-by-1 matrix
%                     MeasureRet(iPoint).Times = [MeasureRet(iPoint).Times;(RunNo-1)*T_Int];
%                     % cumulative sum of frequencies for height estimation
%                     MeasureRet(iPoint).FreqCumSum = sum(MeasureRet(iPoint).Freq);
%                     % array of cured heights
%                     z = handles.cp.Wavelength*T_Int/(2*(handles.cp.n_m-handles.cp.n_L))*MeasureRet(iPoint).FreqCumSum;
%                     MeasureRet(iPoint).Heights = [MeasureRet(iPoint).Heights;z];
                 

                end
            
                
                MeasureRet = MeasureRet_par;
%                 delete(gcp);
%%------- end of parallel computing ----------------
                
            %% Mark the end of exposed curing, i.e., start of dark curing
            ...Compute & Display the exposed and dark curing height which is average across ROI
                         if ((~isempty(RunNo_uvClose))&&(RunNo <= RunNo_uvClose))||(isempty(RunNo_uvClose))
                            for iPoint = 1:nPOI
                                MeasureRet(iPoint).zExposed = MeasureRet(iPoint).Heights(end);
                            end
                            
                            % average using robustfit to remove outliers
                            if length(MeasureRet) < 3
                                zExposedNorminal = mean([MeasureRet.zExposed]);
                            else
%                                 zExposedNorminal_fit= robustfit(1:1:length([MeasureRet.zExposed]),[MeasureRet.zExposed]);
                                zExposedNorminal_fit= robustfit(ones(length([MeasureRet.zExposed]),1),[MeasureRet.zExposed]);
                                zExposedNorminal = max(0,zExposedNorminal_fit(1));
                            end
                            set(handles.ed_ExposedCuredHeight,'String',zExposedNorminal);
                            
                         elseif RunNo > RunNo_uvClose
                            for iPoint = 1:nPOI
                                MeasureRet(iPoint).zDark = MeasureRet(iPoint).Heights(end)- MeasureRet(iPoint).zExposed;
                            end
                            
                            % average using robustfit to remove outliers
                            if length(MeasureRet) < 3
                                zDarkNorminal = mean([MeasureRet.zDark]);
                            else
%                                 zDarkNorminal_fit= robustfit(1:1:length([MeasureRet.zDark]),[MeasureRet.zDark]);                                
                                zDarkNorminal_fit= robustfit(ones(length([MeasureRet.zDark]),1),[MeasureRet.zDark]);
                                zDarkNorminal = max(0,zDarkNorminal_fit(1));
                            end
                            set(handles.ed_DarkCuredHeight,'String',zDarkNorminal);
                        end
                
                % display average measurement results for ROI
                z_All = [MeasureRet.Heights]; % (RunNO+1)-by-iPoint matrix
%                 zNorminal = mean(z_All(end,:)); 
                % average using robustfit to remove outliers
                if length(MeasureRet) < 3
                    zNorminal = mean(z_All(end,:));
                    meanPhase2Pi = mean([MeasureRet.Phase2Pi]);
                else
%                     z_Mean_fit = robustfit(1:1:length(z_All(end,:)),z_All(end,:));
                    zNorminal_fit = robustfit(ones(length(z_All(end,:)),1),z_All(end,:));
                    zNorminal = max(0,zNorminal_fit(1));
%                     meanPhase2Pi_fit = robustfit(1:1:length([MeasureRet.Phase2Pi]),[MeasureRet.Phase2Pi]);
                    meanPhase2Pi_fit = robustfit(ones(length([MeasureRet.Phase2Pi]),1),[MeasureRet.Phase2Pi]);
                    meanPhase2Pi = max(0,meanPhase2Pi_fit(1));
                end
                set(handles.ed_ICM_MeasuredHeight, 'String', sprintf('%.3f',zNorminal));
                set(handles.ed_Phase2Pi,'String',meanPhase2Pi);
            end % end of the data analysis
      end
        guidata(hObject, handles);
        pause(1/handles.cp.FPS);
end

%% Save and report Offline ICM&M Result if measured
if exist('zNorminal','var') % only do this section when measurement was performed
        set(handles.st_InterferogramStatusBar,'String','Saving results of Offline ICM&M measurement and analysis');
        
        for i = 1:length(MeasureRet)
            MeasureRet(i).FittedCoeffs = array2table(MeasureRet(i).FittedCoeffs, 'VariableNames',...
                {'status', 'rsquare', 'I0', 'I1','freqW', 'freq', 'movingHorizon', 'halfLife'});
        end
        MeasParameters = handles.cp;
        save(strcat(handles.cp.ResultFolder,strcat('\Offline_ICM_',datestr(now,'yyyymmdd_HHMMSS'),sprintf('_H%03d_W%03d',g_POI(1, 1),g_POI(2, 1)),'.mat')),...
            'MeasureRet','RunNo_uvClose','meanPhase2Pi',...
            'zNorminal','zExposedNorminal','zDarkNorminal','MeasBeginFrame','MeasParameters',...
            'NumThresholdPixels_Array','NumDarkPixels_Array');% added this line on Aug-06-2016
        
        reportMeasureRet_Offline(MeasureRet,RunNo_uvClose); % save all points Measurements plots
        
        set(handles.st_InterferogramStatusBar,'String','Offline ICM&M results are saved already');
end

%% if "Stop" is pressed, toggle button turns off
if handles.CFrameInd ~= 1 % stop playing so need to go back to current frame instead of upcoming frame
    handles.CFrameInd = handles.CFrameInd-1;
end
if handles.CFrameInd == handles.AVInFrame % stop automatically at last frame
    set(handles.tb_PlayStopICMVideo,'Value',0, 'String', 'Play','Enable', 'on');
end
set(handles.tb_PlayStopICMVideo, 'String', 'Play','Enable', 'on'); % Stop button pushed, reset to "Play"
set(handles.pb_ICMVideo1stFrame, 'Enable', 'on');
set(handles.pb_ICMVideoPrevFrame, 'Enable', 'on');
set(handles.pb_ICMVideoNextFrame, 'Enable', 'on');
set(handles.pb_ICMVideoEndFrame, 'Enable', 'on');
guidata(hObject, handles);
