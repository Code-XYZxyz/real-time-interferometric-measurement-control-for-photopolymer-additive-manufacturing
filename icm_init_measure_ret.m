function MeasStruct = icm_init_measure_ret(cp)

%% Initialize the structure array of measurement result for the chosen POI
    ...Will be used in the "Real_Time_ICM_processMeasureTimer.m"

% data structure for each point measurement with all the fields above
MeasStruct = struct(...
    'PixelHeightWidth',[],... % point coordination [height;width]
    'dataX',[],... % point data: array of time (s) in curve fitting
    'rawY',[],... % point data: time-series of pixel grayscale
    'dataY',[],...
    'fitY',[],... % fitted data using curve fitting
    'FittedCoeffs',[],... % Curve fitting coefficients at initial run per point
    'CureFlags', struct('CureFlag_RunNo',0,'CureFlag_FrameIdx',0),...% Flag the curing window, i.e., mark the beginning of curing
...RunNo and frame index that mark the begin of curing
    'Idx_FailFit',zeros(1,1),... % Flag the runs of failed curve fitting, which has low R-square and may yield frequency outlier Idx_FailFit = [failed RunNO]
    'Times',[0],...% array of run time per point, will grow to be a RunNo-by-1 matrix
    'Heights',[0],...% height estimation at initial run per point, will grow to be a RunNo-by-1 matrix
    'zExposed',[],...% ICM Measured Cured Heights when UV closes
    'zDark',[],...% ICM Measured Cured Heights after UV closes at the end of acquisition
    'Freq_w',[],...% frequency of Im = I0 + I1*cos(wt+phi),Freq_w: angular frequency in Im=I0+I1*cos(W*t+phi)
    'Freq',[],...%Freq: frequency = Freq_w/2/pi (unit:Hz)
    'Phase2Pi',0,...%  the phase (unit: 2Pi), i.e.,time cumulative sum of frequencies to measure height
    'lastFitRet', icm_init_fit_ret(cp)...
    );

% FittedCoeffs: Curve fitting coefficients at initial run per point
..."fourier1" returns 4 coefficients y=a0+a1*cos(px)+b1*sin(px)=I0+I1*cos(Wt+phi)
    ...note: first set of coefficients is the initial values
    ...The length of coefficients = 4*(RunNo + 1) due to the initial value
...FittedCoeffs = [fitStatus,fitgof.rsquare,... % goodness of fitting
    ...I0 (baseline amplitude),I1 (Oscillation amplitude), freqW, freq (Hz),... % estimated parameters
    ...movingHorizon, halfLife], % fitting parameters

