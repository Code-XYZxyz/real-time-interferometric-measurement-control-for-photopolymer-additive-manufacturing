function ret = icmFit2(dataX, dataY, params, prevRet)
%% Function returns the following parameters:
...ret.fitobject = fitobjects{iFitTrial}; % fitting trial which outstands
    ...ret.fitgof = gofs{iFitTrial};
    ...ret.movingHorizon = movingHorizons(iFitTrial);
    ...ret.halfLife = halfLifes(iFitTrial);
    ...ret.frameIdx = frameIdx_local = length(rawY);
    ...ret.firstValidFoiIdx: curing starts,i.e.,end of threshold
    ...ret.fitStatus: 0:not started; 1:fail; 2:valid
    ...ret.fittype: "fourier 1"
    ...ret.fitoptions
    ...ret.I0: baseline amplitude (DC)
    ...ret.I1: oscillation amplitude (AC)
    ...ret.freqW: angular frequency in Im=I0+I1*cos(W*t+phi)
    ...ret.freq: frequency = W/2/pi (unit:Hz)
    ...ret.time: relative time (s) from beginning of measurement = dataX(end)
    %% Curve fitting parameters
% params.rSquare = handles.GOF_rSquare;
% params.FPS = handles.cp.FPS;
% params.halfLife = handles.HalfLife;
% params.MHL = handles.MovingHorizonL;
% params.MeasPeriod = handles.MeasPeriodSamples;
% params.uvIris = uvIris;
% params.RunNo_uvClose = RunNo_uvClose;
% params.RunNo= RunNo;
% params.frameIdx = frameIdx;
% params.f_max = handles.f_max;
% params.f_diff_max = handles.f_diff_max;
%% Set up fitting options: limits, start points
... DONE in this section:
    % 1. Use rsquare to check fit goodness
% 2. Use previous good parameter as starting point of next fit
... TO-DO in next Section:
    % 1. Adaptive half life and moving horizon if rsquare too low


% fprintf('icmFit2 inside starts\n');    

%% Smooth time sequence of grayscale
dataY = smooth(dataY);
dataLen = length(dataY);
% frameIdx_local = params.frameIdx;
frameIdx_local = dataLen;

% fprintf('smooth done\n');    

%% set critical values to identify threshold and dark period
% Threshold_I1 = 10;
% % 1.5s range to detect dark curing
% darkWindowLen = 44; 
% darkRange = 20;

%%--if want to use different criterions for threshold and dark period, use below
if isempty(params.uvIris)|| (params.uvIris >= 10)
    Threshold_I1 = 10;
    % 1.5s range to detect dark curing
    darkWindowLen = 44; 
    darkRange = 20;
else % updated on 07-22-2016: set empirical frequency range [0.1 Hz,0.5 Hz]for uvIris=5
    Threshold_I1 = 10;
    % 2.5s range to detect dark curing
    darkWindowLen = 44; 
    darkRange = 10;
end

%%---------------------

%% Initialize fitting options
if params.RunNo == 1
    
    % y = a0 + a1*cos(w*x) + b1*sin(w*x),where x is actually time "t"
    % Coefficients = [a0;a1;b1;w];
    opts = fitoptions('fourier1');
    opts.Display = 'Off';
    
    % note: setting reasonable limit is very important!
    % set limits of the coefficients "a0,a1,b1,w"
    % If FPS of video is 30, maximum detectable frequency is 30/2=15, while
    %     wMax = 2 * pi * params.FPS/2;
%     wMax = 2 * pi * min(params.FPS/2,params.f_max);
    aMax = 255 * 1;
    
    %%--- 1st curve fitting, bound frequency to get a good start, starting point is important
    % updated on 05-26-2016: set empirical frequency range [0.4 Hz,1.2 Hz] for uvIris>5
    if isempty(params.uvIris)|| (params.uvIris >= 10)
        opts.Lower = [0 -aMax/2 -aMax/2 0.4*2*pi];
        opts.Upper = [aMax aMax/2 aMax/2 1.2*2*pi];
    else % updated on 07-22-2016: set empirical frequency range [0.1 Hz,0.5 Hz]for uvIris=5
        opts.Lower = [0 -aMax/2 -aMax/2 0.1*2*pi]; 
        opts.Upper = [aMax aMax/2 aMax/2 0.5*2*pi];
    end
    %%---------
    
    opts.MaxIter = 400;
    %     opts.TolFun = 1.0e-6;
    %     opts.TolX = 1.0e-6;
    opts.TolFun = 1.0e-3;
    opts.TolX = 1.0e-5;
    % opts.StartPoint = [0 0 0 0];
    
else % if previous fitting exists, keep using it
    opts = prevRet.fitoptions;
    % note: setting reasonable limit is very important!
    % set limits of the coefficients "a0,a1,b1,w"
    % If FPS of video is 30, maximum detectable frequency is 30/2=15, while
    wMax = 2 * pi *min(params.FPS/2,params.f_max);
        %             wMax = 2 * pi *params.f_max;
    aMax = 255 * 1;
    
    % updated on 07-22-2016
    % Set frequency range [0.1,15 Hz] for curing period before UV closes
    if (prevRet.fitStatus ~= 0) && (isempty(params.RunNo_uvClose) || (params.RunNo <= params.RunNo_uvClose))
        opts.Lower = [0 -aMax/2 -aMax/2 0.1*2*pi];
        opts.Upper = [aMax aMax/2 aMax/2 wMax];
    elseif ~isempty(params.RunNo_uvClose) && (params.RunNo > params.RunNo_uvClose) % dark curing period, fitting frequency could be low to 0       
        opts.Lower = [0 -aMax/2 -aMax/2 0];
        opts.Upper = [aMax aMax/2 aMax/2 wMax];
    end
    
    % use the previous starting point to save time
%     if prevRet.fitStatus == 2 || prevRet.fitStatus == 42 ||(prevRet.fitgof.rsquare >= 0.75 && prevRet.fitgof.rsquare ~= 1) % till 07-22-2016
    if prevRet.fitStatus == 2 || prevRet.fitStatus == 42 || (prevRet.fitStatus ~= 0 && prevRet.fitgof.rsquare >= 0.75 && prevRet.fitgof.rsquare ~= 1) % updated 07-22-2016
        opts.StartPoint = coeffvalues(prevRet.fitobject);
    end
end

% fprintf('fit opts done\n');

%% At threshold period, useless fitting, simply return
if (prevRet.fitStatus == 0) && (std(dataY(prevRet.frameIdx+1:end)) < 5)&&(range(dataY(1:end)) < 20)% previously use 10
    ret.fitobject = [];
    ret.fitgof.rsquare = 1;
    ret.movingHorizon = params.MHL;
    ret.halfLife = params.halfLife;
    ret.frameIdx_dummy = frameIdx_local;
    ret.frameIdx = params.frameIdx;
    ret.firstValidFoiIdx = 0; % NOT curing frame yet
    ret.fitStatus = 0; % NO fitting yet
    ret.fitoptions = opts;
    % retrun DC (Direct Current) values
    ret.I0 = mean(dataY((prevRet.frameIdx+1):end));
    ret.I1 = 0;
    ret.freqW = 0;
    ret.freq = 0;
    ret.time = dataX(end);
    return
end

%% At the dark curing period, when the curve flats out, no need fitting
if (~isempty(params.RunNo_uvClose) && (std(dataY(prevRet.frameIdx+1:end)) < 5)...
        &&(range(dataY(end-darkWindowLen:end)) <= darkRange))...%% updated on 07-22-2016
        || (prevRet.fitStatus == 3)
    ret = prevRet;
    ret.fitobject = [];
    ret.fitgof.rsquare = 1;
    ret.movingHorizon = params.MHL;
    ret.halfLife = params.halfLife;
    ret.frameIdx_dummy = frameIdx_local;
    ret.frameIdx = params.frameIdx;
    ret.fitStatus = 3; % Dark curing, NO fitting
    ret.fitoptions = opts;
    % retrun DC (Direct Current) values
    ret.I0 = mean(dataY((prevRet.frameIdx+1):end));
    ret.I1 = 0;
    ret.freqW = 0;
    ret.freq = 0;
    ret.time = dataX(end);
    return
end

%% Start fitting
% R-square threhold value for fitting goodness
rsqTh = params.rSquare;

%% ------------ use window length of approximately one cycle, and 1/3 - 2/3 halflife
% nMaxTrial = 4;
%
% halfLifes = [1.5*params.MeasPeriod, round(params.halfLife* (1:0.5:1+(nMaxTrial-2)*0.5))];
% movingHorizons = min(dataLen,[2*params.MeasPeriod, params.MHL*ones(1,nMaxTrial-1)]);

nMaxTrial = 2;
halfLifes = params.halfLife* [1,1];
movingHorizons = min(dataLen,round(params.MHL*[1,1.5]));
% halfLifes = round(params.halfLife* (1:0.5:2));
% movingHorizons = min(dataLen,params.MHL*ones(1,nMaxTrial));
% if dark curing period, extend the window length to estimate the lower
% frequency
if ~isempty(params.RunNo_uvClose) && (params.RunNo > params.RunNo_uvClose)
    halfLifes = round(params.halfLife* [2,3]);
    movingHorizons = min(dataLen,round(params.MHL*[2,3]));
    %     halfLifes = round(params.halfLife* (1.5:0.5:2.5));
    %     movingHorizons = min(dataLen,2*params.MHL*ones(1,nMaxTrial));
end

% uvIris level too low, e.g., uvIris=5, curing frequency is 0.2Hz, need
...longer MHL to estimate the frequency accurately
    % if (params.uvIris < 10) && (prevRet.fitStatus ~= 0)
if params.uvIris < 10
    halfLifes = params.halfLife* [1,2];
%     movingHorizons = min(dataLen,round(params.MHL*[2,2.5])); % used before 07-23-16 
    movingHorizons = min(dataLen,round(params.MHL*[2,3])); 
end


% halfLifes = [max(10,params.MeasPeriod), round(params.halfLife* (1:0.5:1+(nMaxTrial-2)*0.5))];
% movingHorizons = min(dataLen,[2*max(10,params.MeasPeriod), params.MHL*ones(1,nMaxTrial-1)]);
%% ---------------  icmFit2_v2.m   ---------------
% %% ---- Determine Half life and Moving Horizon Length for fitting
% %%---- Adaptive half life and moving horizon if rsquare too low
% ... Different halfLife and Moving Horizon adjust methods according to prevous fitting status
% ... Try multiple fittings with various half life and MHL until a good fitting achieved
% ...Note: all the window length cannot exceed the available data length
%
% %%------------- Case 1: at the beginning of curing directly after the end of threshold
% ... use mainly the current Measurement period samples by carefully including some precluding samples
% ... because the threshold samples may mislead the fitting and induce errors
% ... However, if the measurement period is too small, e.g., 5 samples per run of measurement,
% ... the small pool of samples fitting may be dominated by noise and yield irrealistic high curing frequency
% ... Hence, we confine the half life with measurement period or 10 which ever is more to minimize the noise effect
% ... and use at most 2 cycles of measurement samples to reduce the threhold samples effect.
% ... Note that threshold samples are not completely unnecessary, presence of proper amount of threshold samples
% ... will help the fitting process learn better about the pattern and increase the fitting confidence.
% if prevRet.fitStatus == 0
% % nMaxTrial = 3;
% % % halfLifes = [prevRet.halfLife,params.MeasPeriod,round(1.5*params.MeasPeriod)];
% % % movingHorizons = min(dataLen,[prevRet.movingHorizon,params.MeasPeriod,2*params.MeasPeriod]);
% % halfLifes = [params.MeasPeriod,round(1.5*params.MeasPeriod),round(2*params.MeasPeriod)];
% % movingHorizons = min(dataLen,[2*params.MeasPeriod,2*params.MeasPeriod,2*params.MeasPeriod]);
% nMaxTrial = 2;
% halfLifes = [max(10,params.MeasPeriod),max(10,params.MeasPeriod)];
% movingHorizons = min(dataLen,[round(1.5*max(10,params.MeasPeriod)),2*max(10,params.MeasPeriod)]);
% end
%
% %%------------- Case 2: previous fitting is NOT valid or failed
% % Method 1. if previous fitting failed, assuming a decreasing curing freq, need
% ... firstly, try the current measurement period sample (use 10 if smaller than 10 samples in the period)
% ...then,with initial half life (e.g. 10) and MH(e.g 32),increase half life by step of (1/2 initial value)
% ...lastly, one may also try double the inital moving horizon(e.g. 2*32=64)
%
% % if prevRet.fitStatus == 1 || prevRet.fitStatus == 3 || prevRet.fitStatus == 41 || prevRet.fitStatus == 40
% % nMaxTrial = 3;
% % halfLifes = round(params.halfLife* (1:0.5:2));
% % movingHorizons = min(dataLen,params.MHL*ones(1,nMaxTrial));
% % end
%
% if prevRet.fitStatus == 1 || prevRet.fitStatus == 3 || prevRet.fitStatus == 41 || prevRet.fitStatus == 40
% % nMaxTrial_half = 3;
% % nMaxTrial = 2*nMaxTrial_half+1;
% nMaxTrial = 4;
% halfLifes = [max(10,params.MeasPeriod), round(params.halfLife* (1:0.5:1+(nMaxTrial-2)*0.5))];
% movingHorizons = min(dataLen,[2*max(10,params.MeasPeriod), params.MHL*ones(1,nMaxTrial-1)]);
% end
%
% %%------------- Case 2: previous fitting is valid
% % Method 2. depending on previous freq,
% ...Step 1: use previous half life & MHL first, if succeed, move on
% ...Step 2: otherwise, use previous period as reference to decide new half life & MHL
%     ...increase half life by step of portion (e.g 0.1) of previous period
%     ...but make sure at least the whole current measurement period samples are included which is important in real-time practice.
%     ...use 2 times of half life as MHL or 2 times measurement period, whichever is larger.
% ...Step 3: use twice measurement period as half life, and four times measurement period as window length
%         % moving horizon: as half life grows with previous period,
%         % too-large window introduces lots of previous curing oscillations
%         % that are different from the slow tailing and cause huge errors
%         % in frequency estimation, hence limit to four times measurement
%         % period
%
%
% if prevRet.fitStatus == 2 || prevRet.fitStatus == 42
% steps = 0.1:0.1:1;
% nMaxTrial = length(steps)+1;
% prevPeriod = 1 / prevRet.freq;
% halfLifes_adpt = round(prevPeriod * params.FPS * steps);
% halfLifes = [prevRet.halfLife,min(max(max(10,params.MeasPeriod),halfLifes_adpt),max(10,params.MeasPeriod)*2)];
% movingHorizons = [prevRet.movingHorizon,...
%     min(dataLen, min(max(max(10,params.MeasPeriod),halfLifes_adpt),max(10,params.MeasPeriod)*2)* 2)];
% end
%
% % 3. TODO: depending on previuos freq and recent freq slope

%% Curve fitting with the options above:
% fit type, limits, start point, moving horizon, weights (half life)
fitobjects = cell(nMaxTrial, 1);
gofs = cell(nMaxTrial, 1);

fitStatus = 1;
for iFitTrial = (1:nMaxTrial)
    halfLife = halfLifes(iFitTrial);
    movingHorizon = movingHorizons(iFitTrial);
    
    trainLen = movingHorizon;
    trainIdxs = ((dataLen - trainLen + 1) : dataLen)';
    trainX = dataX(trainIdxs);
    trainY = dataY(trainIdxs);
    
    weights = exp((trainIdxs - trainIdxs(end)) / halfLife * log(2));
    opts.Weights = weights;
    
%     fprintf('start fit call. size(x)=%d size(y)=%d\n', size(trainX), size(trainY));
    [fitobject, gof] = fit( trainX, trainY, 'fourier1', opts );
%     fprintf('end fit call\n');
    fitobjects{iFitTrial} = fitobject;
    gofs{iFitTrial} = gof;
    
    if round(gof.rsquare,2) >= rsqTh
        fitStatus = 2;
        break;
        %     elseif prevRet.firstValidFoiIdx == 0
        %         % We have not found a single valid fit yet, just stop for now
        %         break;
    else
        if iFitTrial > 1 && gof.rsquare < gofs{Idx_PrevTrial}.rsquare
            %             % half life change is not good, try a 2nd window length with same half-life
            %             if(prevRet.fitStatus == 1) && (iFitTrial<=ceil(nMaxTrial/2))
            %                 Idx_PrevTrial = iFitTrial;
            %                 iFitTrial = iFitTrial+nMaxTrial_half;
            %             else % half life change is not good, use previous fit
            iFitTrial = Idx_PrevTrial;
            break;
            %             end
        else
            % otherwise keep going
            Idx_PrevTrial = iFitTrial;
        end
    end
end

%% Return values
ret.fitobject = fitobjects{iFitTrial};
ret.fitgof = gofs{iFitTrial};
ret.movingHorizon = movingHorizons(iFitTrial);
ret.halfLife = halfLifes(iFitTrial);

ret.frameIdx_dummy = frameIdx_local;
ret.frameIdx = params.frameIdx;
ret.firstValidFoiIdx = prevRet.firstValidFoiIdx;

if fitStatus ~= 0 && prevRet.firstValidFoiIdx == 0
    ret.firstValidFoiIdx = prevRet.frameIdx+1;
end

% if fitStatus == 2 && prevRet.firstValidFoiIdx == 0
%     ret.firstValidFoiIdx = frameIdx_local;
% end

coeffs = coeffvalues(fitobject);
I0 = coeffs(1); % estimated baseline amplitude (DC)
I1 = sqrt(coeffs(2)^2+coeffs(3)^2); % estimated fringe amplitude (AC)
freqW = coeffs(4); % estimated angular frequency "w"
freq = coeffs(4) / 2 / pi;% estimated frequency "f" (unit:Hertz)

ret.fitStatus = fitStatus;
ret.fitoptions = opts;

ret.I0 = I0;
ret.I1 = I1;
ret.freqW = freqW;
ret.freq = freq;
ret.time = dataX(end);

%% Frequency outliers detection & correction
...In the threshold period, correct the frequency artificially by using previous
% the following is especially useful in low uvIris when threshold
% period is long and fluctuating, to identify the threhold correctly is
% very important in the overall process parameters estimation
% if (prevRet.fitStatus == 0) && ((ret.I1 < 10)||(ret.freq < 0.05))
if (prevRet.fitStatus == 0) && (ret.I1 < Threshold_I1)
    ret = prevRet; % use previous returned fitting
    ret.fitgof.rsquare = 0;
    ret.frameIdx_dummy = frameIdx_local;
    ret.frameIdx = params.frameIdx;
    ret.fitStatus = 0;
    ret.fitoptions = opts;
    
    % retrun DC (Direct Current) values
    ret.I0 = mean(dataY((prevRet.frameIdx+1):end));
    ret.I1 = 0;
    ret.freqW = 0;
    ret.freq = 0;
    ret.time = dataX(end);
end

% the following is especially useful in low uvIris when threshold
% period is long and fluctuating, to identify the threhold correctly is
% very important in the overall process parameters estimation
if (prevRet.fitStatus == 0) && (ret.freq < 0.1)
    ret.fitStatus = 0;
end

%--- in the curing period, if too noisy, use zero frequency directly.
if (prevRet.fitStatus ~= 0) && (ret.I1 < 5)
    ret = prevRet; % use previous returned fitting
    ret.fitgof.rsquare = 0;
    ret.frameIdx_dummy = frameIdx_local;
    ret.frameIdx = params.frameIdx;
    ret.fitStatus = 40; % outlier frequency due to failed fitting
    
    % retrun DC (Direct Current) values
    ret.I0 = mean(dataY((prevRet.frameIdx+1):end));
    ret.I1 = 0;
    ret.freqW = 0;
    ret.freq = 0;
    ret.time = dataX(end);
end


end
