function ret = icm_main_worker(cp)

disp('icm main worker thread starts');

% curve fitting parameters
% global RunNo % number of runs of rolling fit

%==== Measurement result return: rolling fit coeff, online height estimate
% nPOI-by-1 structure array of all points measurement
% Each point structure has fields: 'PixelHeightWidth', 'rawY','fitY','FittedCoeffs'
... 'CureFlags', 'Idx_FailFit', 'Times', 'Heights', 'Freq_w','Freq','FreqCumSum'
% Initialized at callback function "pb_SetROI_Callback.m"
% global MeasureRet 

%==== Control variables and results
% global tExpStart % stop watch timer of exposure
% global tExpDuration % Actual exposure duration since it starts (s)
% global frameIdx_uvClose
% global RunNo_uvClose 
% global ExpTimeTarget % Target exposure time in stopwatch control of exposure time
% global ExpTimeNorminal % Norminal exposure time when decide to stop UV
% global zExposedNorminal % ICM Measured Average Cured Heights across ROI when UV closes

%% init memory map file for gui session
clearFile = 0;
mmf = icm_init_mem_file(cp, clearFile);

mmfReadIdx = 1;
mmfMaxNumCache = size(mmf.Data, 1);
fprintf('mmfMaxNumCache=%d\n', mmfMaxNumCache);

%% init

RunNo = 0;
RunNo_uvClose = [];

% nPOI-by-1 structure array of all points measurement  
nPOI = size(cp.POI, 2);
MeasStruct = icm_init_measure_ret(cp);
MeasureRet = repmat(MeasStruct,nPOI,1);  
for iPoint = 1:nPOI
    % single point identification by coordinations (height, width)
    MeasureRet(iPoint).PixelHeightWidth = cp.POI(:, iPoint);
end

lastReadFrameIdx = 0;
nCache = 0;
statPts = zeros(2, 5000);

% tell GUI that worker is ready now
mmf.Data(1).status(3) = 1;

%% start parsing
while mmf.Data(1).status(1) == 1 || lastReadFrameIdx < mmf.Data(1).status(2)
    
    if lastReadFrameIdx >= mmf.Data(1).status(2)
        % keep iterating until there's new frame coming
        continue
    end
    
    % If there's new frame, read it
    cacheData = mmf.Data(mmfReadIdx);
    nCache = nCache + 1;
    frameIdx  = cacheData.frameIdx;
    frameTime = cacheData.frameTime;
    frame     = cacheData.frame;
    snapTic   = cacheData.snapTic;
    uvIris    = cacheData.uvIris;
    uvStatus  = cacheData.uvStatus;
    
    % fixme: preallocate
    cacheDataVec(nCache) = cacheData;
    frameTimeVec(nCache) = frameTime;
    
    if frameIdx ~= nCache
        warning('frameIdx %d should equal to nCache %d', frameIdx, nCache);
    end

    % update cache
    mmf.Data(mmfReadIdx).frameIdx(1) = 0;
    lastReadFrameIdx = frameIdx;
    
    t = toc(snapTic);
    fprintf('received frame %d from cache %d, t=%.3f\n', frameIdx, mmfReadIdx, t);
    statPts(:, nCache) = [double(frameIdx); t];
    
    mmfReadIdx = mod(mmfReadIdx, mmfMaxNumCache) + 1;

    %% prepare fit data
    for iPoint = 1:nPOI
        % single point identification by coordinations (height, width)
        h = cp.POI(1, iPoint); % height coordinate
        w = cp.POI(2, iPoint); % width cooridnate
        
        %NEW method: median filtering neighboring pixels for single pixel intensity to reduce noise
        % Note: frame is a matrix of width-by-height, different from
        ...offline read frame data which is rotated
%         area = frame((h-3):(h+3), (w-3):(w+3)); % 7X7 filter
        area = frame((h-2):(h+2), (w-2):(w+2)); % 5X5 filter
        dataY_foi = double(median(area(:)));
        MeasureRet(iPoint).rawY = [MeasureRet(iPoint).rawY; double(frame(h, w))];
        MeasureRet(iPoint).dataY = [MeasureRet(iPoint).dataY; dataY_foi];
    end
        
    %% Rolling fit: when sufficient samples (set by user,e.g. 20)
    ... and when every measurement period arrived
        % cp.MeasPeriodSamples: online update the model parameters every "MeasPeriodSamples" new data is
    ... acquired, and meanwhile predict next set of "MeasPeriodSamples" data.
        % cp.MeasPeriodSamples = str2double(get(cp.ed_MeasPeriodSamples,'String'));
    if (frameIdx < cp.SamplesNumB4Measure) || (mod(frameIdx-cp.SamplesNumB4Measure, cp.MeasPeriodSamples) ~= 0)
        continue
    end
    
    RunNo = RunNo + 1; % Run number of rolling fit&prediction
    fprintf('Run rolling fit no. %d\n', RunNo);
        
    dataX = frameTimeVec' - frameTimeVec(1); % time of each foi relative to start of measurement
    
    % Curve fitting parameters
    params.rSquare    = cp.GOF_rSquare;
    params.FPS        = cp.FPS;
    params.halfLife   = cp.HalfLife;
    params.MHL        = cp.MovingHorizonL;
    params.MeasPeriod = cp.MeasPeriodSamples;
    params.f_max      = cp.f_max;
    params.f_diff_max = cp.f_diff_max;
 
    params.uvIris     = uvIris;
    params.RunNo      = RunNo;
    params.frameIdx   = frameIdx;
    % Mark the end of exposed curing, i.e., start of dark curing
    % if  (MeasBeginFrame + frameIdx >= frameIdx_uvClose) && (isempty(RunNo_uvClose))
    %     RunNo_uvClose = RunNo;% flag the frame number in acquired video when UV closes
    % end
    params.RunNo_uvClose = RunNo_uvClose;

    %% point-by-poiont analysis
%     for iPoint = 1:nPOI
    parfor iPoint = 1:nPOI
        % single point identification by coordinations (height, width)
        h = cp.POI(1, iPoint); % height coordinate
        w = cp.POI(2, iPoint); % width cooridnate
        fprintf('fitting point %d [%d %d]\n', iPoint, h, w);
           
        MeasureRet(iPoint).dataX = dataX;
        
        % time series of intensity data
        dataY = MeasureRet(iPoint).dataY; % without missing data imputation
        
        %% Curve Fitting
        % Rolling fit with "fourier1" returns 4 coefficients y=a0+a1*cos(px)+b1*sin(px)
        % fitRollRet = icmFit(trainX, trainY, trainW);
        % save the fitting coefficients, i.e., online estimates of parameters
        prevFitRet = MeasureRet(iPoint).lastFitRet;
%         fprintf('dataX %d dataY %d\n', size(dataX), size(dataY));
        fitRollRet = icmFit2(dataX, dataY, params, prevFitRet);
        MeasureRet(iPoint).lastFitRet = fitRollRet;
        
        % save the fitting coefficients, i.e., online estimates of parameters
        coeffs = [fitRollRet.fitStatus,fitRollRet.fitgof.rsquare,...
            fitRollRet.I0, fitRollRet.I1, fitRollRet.freqW, fitRollRet.freq,...
            fitRollRet.movingHorizon, fitRollRet.halfLife];
        
        MeasureRet(iPoint).FittedCoeffs = [MeasureRet(iPoint).FittedCoeffs; coeffs];
        MeasureRet(iPoint).Freq_w = MeasureRet(iPoint).FittedCoeffs(:,5);
        MeasureRet(iPoint).Freq = MeasureRet(iPoint).FittedCoeffs(:,6);
        
        if isfield(fitRollRet,'fitobject')
            if (fitRollRet.fitStatus ~=0 && fitRollRet.fitStatus ~=3 && fitRollRet.fitStatus ~=40)...
                    ||(fitRollRet.fitStatus ==0 && fitRollRet.freq~=0 ) % small freq (<0.1Hz)in threshold acceptable
                newFitY = feval(fitRollRet.fitobject, dataX( (prevFitRet.frameIdx+1):fitRollRet.frameIdx));
                %                             newFitY = feval(fitRollRet.fitobject, dataX( end-cp.MeasPeriodSamples+1:end));
                MeasureRet(iPoint).fitY = [MeasureRet(iPoint).fitY; newFitY];
            else
                newFitY = fitRollRet.I0*ones(fitRollRet.frameIdx-prevFitRet.frameIdx,1);
                MeasureRet(iPoint).fitY = [MeasureRet(iPoint).fitY; newFitY];
            end
        else
            newFitY = fitRollRet.I0*ones(fitRollRet.frameIdx-prevFitRet.frameIdx,1);
            MeasureRet(iPoint).fitY = [MeasureRet(iPoint).fitY; newFitY];
        end
        
        % Mark the start of curing, i.e., the end of threshold
        if  (MeasureRet(iPoint).CureFlags.CureFlag_RunNo==0) && (fitRollRet.firstValidFoiIdx~=0)
            MeasureRet(iPoint).CureFlags.CureFlag_RunNo = RunNo;
            MeasureRet(iPoint).CureFlags.CureFlag_FrameIdx = fitRollRet.firstValidFoiIdx;
        end
        % Flag the runs of failed curve fitting, which has low R-square and may yield frequency outlier
        if fitRollRet.fitStatus ~= 2
            MeasureRet(iPoint).Idx_FailFit = [MeasureRet(iPoint).Idx_FailFit; RunNo];
        end
            
        %% Height Estimation: growth by integration
        T_Int = dataX(end)- prevFitRet.time;
        % array of measurement time(s) per point, RunNo-by-1 matrix
        MeasureRet(iPoint).Times = [MeasureRet(iPoint).Times;dataX(end)];
        
        % phase(unit: 2Pi: time cumulative sum of frequency-by-time for height estimation
        %                    %--- method 1: mixed use of trapzoidal and local value
        %                    if T_Int < 1.5*cp.MeasPeriodSamples/cp.FPS
        %                        MeasureRet(iPoint).Phase2Pi = MeasureRet(iPoint).Phase2Pi+T_Int*fitRollRet.freq;
        %                    else  % if too long interval, use average freq
        %                        MeasureRet(iPoint).Phase2Pi = MeasureRet(iPoint).Phase2Pi+T_Int*(fitRollRet.freq+prevFitRet.freq)/2;
        %                    end
        %--- method 2: always use trapzoidal rule, i.e., midpoint
        MeasureRet(iPoint).Phase2Pi = MeasureRet(iPoint).Phase2Pi+T_Int*(fitRollRet.freq+prevFitRet.freq)/2;
        %                    %--- method 3: always use local frequency only
        %                    MeasureRet(iPoint).Phase2Pi = MeasureRet(iPoint).Phase2Pi+T_Int*fitRollRet.freq;
        
        % array of cured heights
%         % Before 08/24/2016: use constant refractive index
%         z = cp.Wavelength/(2*(cp.n_m-cp.n_L))*MeasureRet(iPoint).Phase2Pi;
        
        % Created on 08/24/2016: use evolving refractive index
        n_m_evolve = 0.00041*(MeasureRet(iPoint).Phase2Pi)+1.49191;% 08/24/2016:calculate evolving refractive index with the model in thesis
        z = cp.Wavelength/(2*(n_m_evolve-cp.n_L))*MeasureRet(iPoint).Phase2Pi; % 08/24/2016: evolving refractive index
        % end of 08/24 updates
        
        MeasureRet(iPoint).Heights = [MeasureRet(iPoint).Heights;z];
        
            
    end % end of point-by-poiont analysis
    
    %% Evaluate the average height of a line profile
    heightsALL =[MeasureRet.Heights];
    heights = heightsALL(end,:);
%     % Method 1: quantile average
%     ql = quantile(heights, 0.25);
%     qh = quantile(heights, 0.75);
%     filteredHeights = heights(heights >= ql & heights <= qh);
%     zNorminal = mean(filteredHeights); % average height evaluated for the line

    % Method 2: robustfit
    if length(MeasureRet) < 3
        zNorminal = mean(heights);
        meanPhase2Pi = mean([MeasureRet.Phase2Pi]);
    else
        zNorminal_fit= robustfit(ones(length(heights),1),heights);
        zNorminal = max(0, zNorminal_fit(1));
        meanPhase2Pi_fit = robustfit(ones(length([MeasureRet.Phase2Pi]),1),[MeasureRet.Phase2Pi]);
        meanPhase2Pi = max(0, meanPhase2Pi_fit(1));
    end

    
% set(cp.ed_ICM_MeasuredHeight, 'String', sprintf('%.3f',z_Mean));
% set(cp.ed_Phase2Pi,'String',meanPhase2Pi);
    
    
    mmf.Data(1).avgTotalHeight(1) = zNorminal;
    mmf.Data(1).avgTotalPhase(1) = meanPhase2Pi;

    %% control after receiving new frame & Exposed Curing
    % Real-time control:stopwatch for target exposure time or target height
    % Turn off UV immediately when measured time or height hits target

    if uvStatus == 1
        % Simple  control stopwatch time control
        if cp.targetMode == 2 && ((cp.targetExpTime - frameTime) <= cp.MeasPeriodSamples*0.03)
            sprintf('frameTime %.2f exceeds targetExpTime %.2f, sending msg to shutdown uv', frameTime, cp.targetExpTime)
            mmf.Data(1).status(3) = 2;   
            R.targetExpTime = cp.targetExpTime;  % setpoint of epxosure time
        
            % calculate the total exposure time
%             tExpDuration = dataX(end);  % before 8/26
            tExpDuration_Ideal = frameTime;  % 08/26/2016
%             % to-do: display exposure time in GUI
%             set(cp.ed_ExposureTime, 'string', tExpDuration);
            
            %%---- Mark the end of exposed curing 
            frameIdx_uvClose_Ideal = frameIdx; % flag the frame number in acquired video when UV closes
            RunNo_uvClose = RunNo;
            if exist('MeasureRet','var') && ~isempty(MeasureRet)
                for iPoint = 1:length(MeasureRet)
                    MeasureRet(iPoint).zExposed = MeasureRet(iPoint).Heights(end);
                end
                % average using robustfit to remove outliers
                if length(MeasureRet) < 3
                    zExposedNorminal = mean([MeasureRet.zExposed]);
                else
                    zExposedNorminal_fit= robustfit(ones(length([MeasureRet.zExposed]),1),[MeasureRet.zExposed]);
                    zExposedNorminal = max(0, zExposedNorminal_fit(1));
                end
%                 % to-do: display the exposed cured height in GUI
%                 set(cp.ed_ExposedCuredHeight,'String',zExposedNorminal); 
            end

        end
        
        % Measurement feedback control
        if cp.targetMode == 1 && (cp.targetCuredHeight*0.9-zNorminal <= 0.5)
            sprintf('measured height %.2f exceeds targetHeight %.2f, sending msg to shutdown uv', zNorminal, cp.targetCuredHeight)
            mmf.Data(1).status(3) = 2;
            R.targetCuredHeight = cp.targetCuredHeight; % Save above already: setpoint of cured height
            
            % calculate the total exposure time
%             tExpDuration = dataX(end);  % before 8/26
            tExpDuration_Ideal = frameTime;  % 08/26/2016
%                 % to-do: display exposure time in GUI
%                 set(cp.ed_ExposureTime, 'string', tExpDuration);

            %%---- Mark the end of exposed curing 
            frameIdx_uvClose_Ideal = frameIdx; % flag the frame number in acquired video when UV closes
            RunNo_uvClose = RunNo;

            if exist('MeasureRet','var') && ~isempty(MeasureRet)
                for iPoint = 1:length(MeasureRet)
                    MeasureRet(iPoint).zExposed = MeasureRet(iPoint).Heights(end);
                end
                % average using robustfit to remove outliers
                if length(MeasureRet) < 3
                    zExposedNorminal = mean([MeasureRet.zExposed]);
                else
                    zExposedNorminal_fit= robustfit(ones(length([MeasureRet.zExposed]),1),[MeasureRet.zExposed]);
                    zExposedNorminal = max(0, zExposedNorminal_fit(1));
                end
%                     % to-do: display the exposed cured height in GUI
%                     set(cp.ed_ExposedCuredHeight,'String',zExposedNorminal); 
            end
        end            
       
    end
    %%------------ End of control --------------------------
    
    %% ---------------- Dark curing -------------------%%
    if uvStatus == 2 % display dark curing height
    for iPoint = 1:length(MeasureRet)
       MeasureRet(iPoint).zDark = MeasureRet(iPoint).Heights(end)- MeasureRet(iPoint).zExposed;
    end
    
    % average using robustfit to remove outliers
    if length(MeasureRet) < 3
        zDarkNorminal = mean([MeasureRet.zDark]);
    else
        zDarkNorminal_fit= robustfit(ones(length([MeasureRet.zDark]),1),[MeasureRet.zDark]);
        zDarkNorminal = max(0, zDarkNorminal_fit(1));
    end
%     % to-do: display the dark cured height in GUI
%     set(cp.ed_DarkCuredHeight,'String',zDarkNorminal);
    end
    
   
end

%% Calculates latence time between acquiring and analyzing a frame
statPts = statPts(:, 1:nCache);

frames = statPts(1,:);
meanDelay = mean(statPts(2,:));

R.rtFramesLatence = statPts;
R.rtProcessedFrames = length(frames);
R.rtMeanDelayAcqAnl = meanDelay;

ret = MeasureRet;
fprintf('icm main worker thread stops, processed %d frames, average delay %.3f\n', length(frames), meanDelay);

clear mmf;

% %% Evaluate the average height of a line profile
% heightsALL =[MeasureRet.Heights];
% heights = heightsALL(end,:);
% % Method 1: quantile average
% ql = quantile(heights, 0.25);
% qh = quantile(heights, 0.75);
% filteredHeights = heights(heights >= ql & heights <= qh);
% zNorminal = mean(filteredHeights); % average height evaluated for the line
% 
% % % Method 2: robustfit
% % if length(MeasureRet) < 3
% %     zNorminal = mean(heights);
% % else
% %     zNorminal_fit= robustfit(ones(length(heights),1),heights);
% %     zNorminal = max(0, zNorminal_fit(1));
% % end


%% Saving results
if cp.isRT == 1 % Real-time
    
    fn = strcat(cp.ResultFolder,...
        strcat('\RT_ECPL_ICM_',datestr(now,'yyyymmdd_HHMMSS'),...
        sprintf('_H%03d_W%03d', cp.POI(1, 1), cp.POI(2, 1)),'.mat'));
    fprintf('saving ret to file %s', fn);
    R.MeasureRet = MeasureRet;
    
%     R.targetExpTime = cp.targetExpTime;  % Saved above already:setpoint of epxosure time
   
    R.frameIdx_uvClose_Ideal = frameIdx_uvClose_Ideal;
    uvStatusVec = [cacheDataVec.uvStatus];
    uvStatusVec_diff = diff(uvStatusVec);
    frameIdx_uvClose = find(uvStatusVec_diff~=0)+1;
    R.frameIdx_uvClose = frameIdx_uvClose;
    R.tExpDuration_Ideal = tExpDuration_Ideal; 
    R.tExpDuration = frameTimeVec(frameIdx_uvClose);
    R.RunNo_uvClose = RunNo_uvClose;
    
%     R.targetCuredHeight = cp.targetCuredHeight; % Save above already: setpoint of cured height
    R.zExposedNorminal = zExposedNorminal; % average exposed height evaluated for the line
    R.zDarkNorminal = zDarkNorminal; % average dark height evaluated for the line
    R.zNorminal = zNorminal; % average height evaluated for the line
    R.cp = cp;
    R.dp = cacheDataVec;
    R.imageTime = [cacheDataVec.frameTime]'; % image time which is required in offline ICM
    save(fn, '-struct', 'R');
end

%% Report results
reportMeasureRet_RT(MeasureRet,RunNo_uvClose);
% reportMeasureRet(MeasureRet(1));

%     if cp.isRT == 1 % Real-time
% 
%         save(strcat(cp.ResultFolder,strcat('\RT_ECPL_ICM_',datestr(now,'yyyymmdd_HHMMSS'),sprintf('_H%03d_W%03d', cp.POI(1, 1), cp.POI(2, 1)),'.mat')),...
%             'MeasureRet','uvIris','tExpDuration','ExpTimeTarget','ExpTimeNorminal',...
%             'RunNo_uvClose','zTarget','zExposedNorminal','zDarkNorminal',...
%             'cp.dp.FrameTimeAbs','cp.dp.frameIdx','MeasParameters'); 
%         reportMeasureRet_RT(MeasureRet,RunNo_uvClose); % save all points Measurements plots
%         reportMeasureRet(MeasureRet(1));
% 
%     %% Offline report    
%     else 
%         for i = 1:length(MeasureRet)
%             MeasureRet(i).FittedCoeffs = array2table(MeasureRet(i).FittedCoeffs, 'VariableNames',...
%                 {'status', 'rsquare', 'I0', 'I1','freqW', 'freq', 'movingHorizon', 'halfLife'});
%         end
%         save(strcat(cp.ResultFolder,strcat('\Offline_ICM_',datestr(now,'yyyymmdd_HHMMSS'),sprintf('_H%03d_W%03d',cp.POI(1, 1),cp.POI(2, 1)),'.mat')),...
%             'MeasureRet','RunNo_uvClose','zExposedNorminal','zDarkNorminal','MeasParameters');
%         
%         % Report results
%         reportMeasureRet_Offline(MeasureRet,RunNo_uvClose); % save all points Measurements plots
% %         reportMeasureRet(MeasureRet(1));
% 
%     end

%% Real-time control for single point only: stopwatch for target exposure time
% ...i.e., turn off UV immediately when measured time hits target
% ...note: the code could be easily adapted to multipoint control
% if (cp.dp.uvStatus == 1)&& get(cp.rb_TargetExpTime, 'Value') == 1 % Stopwatch of Exposure time
%     ExpTimeTarget = str2double(get(cp.ed_TargetExpTime,'String'));
%     ExpTimeNorminal = toc(tExpStart);
% 
%     % when UV light is on and target time reached,stop UV
%     if (ExpTimeTarget-ExpTimeNorminal) <= 0.03  
%          % --- Close UV Shutter
% %             UVShutterClose(uv);
% 
%             % calculate the total exposure time
%             tExpDuration = toc(tExpStart); 
%             % display exposure time
%             set(cp.ed_ExposureTime, 'string', tExpDuration);
% 
%             % flag UV light is off again after being "on", not "0" so that acquisition and measurement could go on to capture dark curing
%             cp.dp.uvStatus = 2; 
% 
%             %%---- Mark the end of exposed curing 
%             frameIdx_uvClose = cp.dp.frameIdx; % flag the frame number in acquired video when UV closes
%             RunNo_uvClose = RunNo;
%             if exist('MeasureRet') && ~isempty(MeasureRet)
%                 for iPoint = 1:length(MeasureRet)
%                     MeasureRet(iPoint).zExposed = MeasureRet(iPoint).Heights(end);
%                 end
%                 % average using robustfit to remove outliers
%                 if length(MeasureRet) < 3
%                     zExposedNorminal = mean([MeasureRet.zExposed]);
%                 else
%                     zExposedNorminal_fit= robustfit(1:1:length([MeasureRet.zExposed]),[MeasureRet.zExposed]);
%                     zExposedNorminal = max(0, zExposedNorminal_fit(1));
%                 end
%                 set(cp.ed_ExposedCuredHeight,'String',zExposedNorminal);
%             end
% 
%            % --- Disconnect UV
% %             UVDisc(uv);
% 
%             % Change button string
%             set(cp.pb_OpenCloseUV,'String','Open UV Light');
% 
%             % update the Interferogram Status Bar for info
%         set(cp.st_InterferogramStatusBar,'String','ON Target TIME. Controller closed UV shutter. Stop measurement when ready.');
%     end  
% end

% % append to data structure
% % Note: difference between frameIdx and cp.dp.frameIdx
% ... is that frameIdx marks the dataset (frames) with ROI to be analyzed
% v3d(:,:,frameIdx) = frame;



% %% Real-time display measurement results for center point
% % display average measurement results for ROI
% z_All = [MeasureRet.Heights]; % (RunNO+1)-by-iPoint matrix
% 
% % average using robustfit to remove outliers
% if length(MeasureRet) < 3
%     z_Mean = mean(z_All(end,:));
%     meanPhase2Pi = mean([MeasureRet.Phase2Pi]);
% else
%     z_Mean_fit = robustfit(1:1:length(z_All(end,:)),z_All(end,:));
%     z_Mean = max(0,z_Mean_fit(1));
%     meanPhase2Pi_fit = robustfit(1:1:length([MeasureRet.Phase2Pi]),[MeasureRet.Phase2Pi]);
%     meanPhase2Pi = max(0, meanPhase2Pi_fit(1));
% end
% set(cp.ed_ICM_MeasuredHeight, 'String', sprintf('%.3f',z_Mean));
% set(cp.ed_Phase2Pi,'String',meanPhase2Pi);
% 
% if cp.dp.uvStatus == 2 % display dark curing height
%     for iPoint = 1:length(MeasureRet)
%        MeasureRet(iPoint).zDark = MeasureRet(iPoint).Heights(end)- MeasureRet(iPoint).zExposed;
%     end
%     
%     % average using robustfit to remove outliers
%     if length(MeasureRet) < 3
%         zDarkNorminal = mean([MeasureRet.zDark]);
%     else
%         zDarkNorminal_fit= robustfit(1:1:length([MeasureRet.zDark]),[MeasureRet.zDark]);
%         zDarkNorminal = max(0, zDarkNorminal_fit(1));
%     end
%     set(cp.ed_DarkCuredHeight,'String',zDarkNorminal);
% end

% %% Real-time control for single point only: stopwatch for target cured height
% ...i.e., turn off UV immediately when measured height hits target
% ...note: the code could be easily adapted to multipoint control
% if (cp.dp.uvStatus == 1)&& get(cp.rb_TargetCuredHeight, 'Value') == 1 % Stopwatch of Cured height
%     zTarget = str2double(get(cp.ed_TargetCuredHeight,'String'));
%     z_All = [MeasureRet.Heights]; % (RunNO+1)-by-iPoint matrix
%     % ICM Measured Average Cured Heights across ROI, used to decide the time to close UV
%     % average using robustfit to remove outliers
%     if length(MeasureRet) < 3
%         zExposedNorminal = mean(z_All(end,:));
%     else
%         zExposedNorminal_fit= robustfit(1:1:length(z_All(end,:)),z_All(end,:));
%         zExposedNorminal = max(0, zExposedNorminal_fit(1));
%     end
%     set(cp.ed_ExposedCuredHeight,'String',zExposedNorminal);
%     
%     % calculate the total exposure time
%     tExpDuration = toc(tExpStart); 
%     % display exposure time
%     set(cp.ed_ExposureTime, 'string', tExpDuration);
% 
%     % when UV light is on and 6/7target height reached,stop UV
%     ...(because dark curing contribute about 1/7 final height)
%     if round(zExposedNorminal) >= zTarget 
%     
%          %-------- Close UV Shutter
% %             UVShutterClose(uv);
% 
%             % calculate the total exposure time
%             tExpDuration = toc(tExpStart); 
%             % display exposure time
%             set(cp.ed_ExposureTime, 'string', tExpDuration);
% 
%             % flag UV light is off again after being "on", not "0" so that acquisition and measurement could go on to capture dark curing
%             cp.dp.uvStatus = 2; 
% 
% 
%             %%---- Mark the end of exposed curing 
%             frameIdx_uvClose = cp.dp.frameIdx; % flag the frame number in acquired video when UV closes
%             RunNo_uvClose = RunNo;
%             for iPoint = 1:nPOI
%                 MeasureRet(iPoint).zExposed = MeasureRet(iPoint).Heights(end);
%             end
% %             zExposedNorminal = mean([MeasureRet.zExposed]);
% %             set(cp.ed_ExposedCuredHeight,'String',zExposedNorminal);
% 
%  
%             %-------- Disconnect UV
% %             UVDisc(uv);
% 
%             % Change button string
%             set(cp.pb_OpenCloseUV,'String','Open UV Light');
% 
%             % update the Interferogram Status Bar for info
%         set(cp.st_InterferogramStatusBar,'String','ON Target HEIGHT. Controller closed UV shutter. Stop measurement when ready.');
%     end    
% end

