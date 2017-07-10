function reportMeasureRet_RT(MeasureRet,RunNo_uvClose)
%% This function
...If use individually, make sure to load required mat files as below:
    ...1)ECPL_***.mat
    ...2)RT_Measurements_***.mat
...1) Reports real-time measurment results
...Note: "mr" is a NumPoints-by-1 structure of measurement return
...2) Calculates the dark curing frequency and contribution to the height estimation
if ~isstruct(MeasureRet)
    return
end

% if structure fields are not empty
if ~istable(MeasureRet(1).FittedCoeffs) 
    for i = 1:length(MeasureRet)
            MeasureRet(i).FittedCoeffs = array2table(MeasureRet(i).FittedCoeffs, 'VariableNames',...
                {'status', 'rsquare', 'I0', 'I1','freqW', 'freq', 'movingHorizon', 'halfLife'});
    end
end

%% result folder
% make a folder for figures storage
if ~exist('Result', 'dir')
mkdir(pwd,'Result');
end
ResultFolder = fullfile(pwd,'Result');
SaveFolder = fullfile(ResultFolder,strcat('RTMeasRetPlots_',datestr(now,'yyyymmdd_HHMMSS')));
if ~exist(SaveFolder, 'dir')
% mkdir(ResultFolder,'RTMeasRetPlots');
mkdir(SaveFolder);
end

%% data
TimesAll=[MeasureRet.Times];
dataX_All = [MeasureRet.dataX];
rawY_All = [MeasureRet.rawY]; % raw pixel grayscale
dataY_All = [MeasureRet.dataY];% filtered Y with median image filter
fitY_All = [MeasureRet.fitY];
HeightsAll=[MeasureRet.Heights];
endHeight = HeightsAll(end,:);
zExposedAll = [MeasureRet.zExposed];
zDarkAll = [MeasureRet.zDark];
Phase2PiAll = [MeasureRet.Phase2Pi];
FreqAll=[MeasureRet.Freq];

%% Grayscales raw data vs. fitted data for All Points: RunNo-by-NumPoints matrix
% save figures of raw and fitted grayscale data
for iPoint = 1:length(MeasureRet)
    figY = figure();
    set(figY, 'visible', 'off');
    PixelH = MeasureRet(iPoint).PixelHeightWidth(1);
    PixelW = MeasureRet(iPoint).PixelHeightWidth(2);
    dataX = dataX_All(:,iPoint);
    rawY = rawY_All(:,iPoint);
    dataY = dataY_All(:,iPoint);
    fitY = fitY_All(:,iPoint);
    % saving raw, filtered, fitted data
    plot(dataX(1:min(length(dataX),length(rawY))), rawY(1:min(length(dataX),length(rawY))),'c-',...
        MeasureRet(iPoint).dataX(1:min(length(MeasureRet(iPoint).dataX),length(MeasureRet(iPoint).dataY))),...
    MeasureRet(iPoint).dataY(1:min(length(MeasureRet(iPoint).dataX),length(MeasureRet(iPoint).dataY))),'b.',...
        dataX(1:min(length(dataX),length(fitY))),fitY(1:min(length(dataX),length(fitY))),'r-');
    legend('Raw Data','Filtered Data for Curve Fitting','Fitted Data');
    
%     % saving only filtered and fitted data 
%     plot(MeasureRet(iPoint).dataX(1:min(length(MeasureRet(iPoint).dataX),length(MeasureRet(iPoint).dataY))),...
%     MeasureRet(iPoint).dataY(1:min(length(MeasureRet(iPoint).dataX),length(MeasureRet(iPoint).dataY))),'b.',...
%         dataX(1:min(length(dataX),length(fitY))),fitY(1:min(length(dataX),length(fitY))),'r-');
%     legend('Filtered Data for Curve Fitting','Fitted Data');
    
    title({'Real-time ICM&M Time Series of Grayscale',...
       sprintf('Pixel (Height,Width):(%03d, %03d)', PixelH, PixelW),...
      sprintf('Estimated Total Phase Angle: %0.4f Cycles; Cured Height: %0.3f um', Phase2PiAll(iPoint),endHeight(iPoint))});
    xlabel('Time(s)');ylabel('Grayscale');
    axis([0 max(dataX)+2 0 300]);
    
    pointYfig = sprintf('%s\\RT_Grayscale_Px_%03d_%03d.png', SaveFolder, PixelH, PixelW);
    saveas(figY, pointYfig);    
    close(figY);
    
end
% %% save figures of Heights for All Points: RunNo-by-NumPoints matrix
% % save figures
% for iPoint = 1:length(MeasureRet)
%     figZ = figure();
%     set(figZ, 'visible', 'off');
%     PixelH = MeasureRet(iPoint).PixelHeightWidth(1);
%     PixelW = MeasureRet(iPoint).PixelHeightWidth(2);
%     t = TimesAll(:,iPoint);
%     z = HeightsAll(:,iPoint);
%     z_uvClose = MeasureRet(iPoint).zExposed;
%     t_uvClose = MeasureRet(iPoint).Times(RunNo_uvClose+1);
%     plot(t, z);
%     hold on
%     plot(t,z_uvClose*ones(length(t),1),'r--');
%     hold on
%     plot(t_uvClose*ones(length(z),1),z,'g--');
%     title({'Real-time ICM&M Estimated Cured Height',...
%        sprintf('Pixel (Height,Width):(%03d, %03d)', PixelH, PixelW) });
%     xlabel('time (s)');ylabel('Cured Height (\mum)');
%     axis([0 max(TimesAll(end,:))+2 0 round(max(HeightsAll(end,:)))+20]);
%     legend(sprintf('Estimated Final Total Cured Height: %0.3f', z(end)),...
%        sprintf('Estimated Exposed Cured Height: %0.3f', z_uvClose) ,...
%        sprintf('Approximate Time when UV Lamp Closed: %0.3f', t_uvClose));
%     pointZfig = sprintf('%s\\RT_ICMM_Height_Px_%03d_%03d.png', SaveFolder, PixelH, PixelW);
%     saveas(figZ, pointZfig);    
%     close(figZ);
%     
% end
% 
% %% save figures of Frequency for All Points: RunNo-by-NumPoints matrix
% % save figures
% for iPoint = 1:length(MeasureRet)
%     figF = figure();
%     set(figF, 'visible', 'off');
%     PixelH = MeasureRet(iPoint).PixelHeightWidth(1);
%     PixelW = MeasureRet(iPoint).PixelHeightWidth(2);
%     t = TimesAll(:,iPoint);
%     freq = FreqAll(:,iPoint);
%     plot(t, [0;freq]);
%     title({'Real-time ICM&M Estimated Instantaneous Frequency (Hz)',...
%        sprintf('Pixel (Height,Width):(%03d, %03d)', PixelH, PixelW) });
%     xlabel('time (s)');ylabel('Instantaneous Frequency (Hz)');
%     axis([0 max(TimesAll(end,:))+2 0 2]);
%     legend(sprintf('Estimated Total Phase Angle: %0.4f Cycles', Phase2PiAll(iPoint)));
%     pointFfig = sprintf('%s\\RT_ICMM_Freq_Px_%03d_%03d.png', SaveFolder, PixelH, PixelW);
%     saveas(figF, pointFfig);    
%     close(figF);
%     
% end

%% plot the estimated height and phase for all the pixels in the ROI
% 2-by-NumPoints matrix, Row 1 Heights, Row 2 Widths
Pixels = [MeasureRet.PixelHeightWidth]; 
HR = Pixels(1,1):5:Pixels(1,end);
WR = Pixels(2,1):5:Pixels(2,end);
lenH = length(HR);
lenW = length(WR);
m = zeros(lenH, lenW);
zExposed = zeros(lenH, lenW);
zDark = zeros(lenH, lenW);
phase2Pi = zeros(lenH, lenW);
for iPoint = 1:length(MeasureRet)
    h = Pixels(1,iPoint);
    w = Pixels(2,iPoint);

    ih = (h - Pixels(1,1))/5 + 1;
    iw = (w - Pixels(2,1))/5 + 1;

    m(ih, iw) = HeightsAll(end,iPoint);
    zExposed(ih, iw) = zExposedAll(iPoint);
    zDark(ih, iw) = zDarkAll(iPoint);
    phase2Pi(ih, iw) = Phase2PiAll(iPoint);
end

% Calculate and save the statistics of measured pixels
if length(MeasureRet)<3
    % Average height
%     endHeight = HeightsAll(end,:);
    meanHeight_ls = mean(HeightsAll(end,:));
    sigmaHeight_ls = std(HeightsAll(end,:)); % rmse from least squares fit
    
    mean_zExposed_ls = mean(zExposedAll);
    sigma_zExposed_ls = std(zExposedAll); % rmse from least squares fit
    
    mean_zDark_ls = mean(zDarkAll);
    sigma_zDark_ls = std(zDarkAll); % rmse from least squares fit
    
    % Average phase
    meanPhase2Pi_ls = mean(Phase2PiAll);
    sigmaPhase2Pi_ls = std(Phase2PiAll); % rmse from least squares fit
    
    save(strcat(SaveFolder,'\MeasStats.mat'), 'endHeight','meanHeight_ls','sigmaHeight_ls',...
    'mean_zExposed_ls','sigma_zExposed_ls',...
    'mean_zDark_ls','sigma_zDark_ls',...
    'meanPhase2Pi_ls','sigmaPhase2Pi_ls');

    % for plot, to be consistent with robustfit results plot, use same terms
    meanHeight = meanHeight_ls;
    sigmaHeight = sigmaHeight_ls;    
    mean_zExposed = mean_zExposed_ls;
    sigma_zExposed = sigma_zExposed_ls;     
    mean_zDark = mean_zDark_ls;
    sigma_zDark = sigma_zDark_ls;     
    meanPhase2Pi = meanPhase2Pi_ls;
    sigmaPhase2Pi = sigmaPhase2Pi_ls;
else
    % Average height using "robustfit" to remove outlier
%     meanHeight_fit = robustfit(1:1:length(HeightsAll(end,:)),HeightsAll(end,:));
    [meanHeight_fit,Height_stats] = robustfit(ones(length(HeightsAll(end,:)),1),HeightsAll(end,:));
    meanHeight = max(0,meanHeight_fit(1)); % robust fit mean
    meanHeight_SE = Height_stats.se(1); % standard error of meanHeight estimates
    sigmaHeight = Height_stats.s; % robust fit sigma
    meanHeight_ls = mean(HeightsAll(end,:)); % ordinary mean
    sigmaHeight_ls = Height_stats.ols_s; % rmse from least squares fit
%     endHeight = HeightsAll(end,:);
    Height_3sigma_Num=endHeight(find(endHeight>=(meanHeight-3*sigmaHeight)& endHeight<=(meanHeight+3*sigmaHeight)));
    Height_3sigma_Percent = length(Height_3sigma_Num) / length(endHeight);
    Height_2sigma_Num=endHeight(find(endHeight>=(meanHeight-2*sigmaHeight)& endHeight<=(meanHeight+2*sigmaHeight)));
    Height_2sigma_Percent = length(Height_2sigma_Num) / length(endHeight);
    Height_1sigma_Num=endHeight(find(endHeight>=(meanHeight-sigmaHeight)& endHeight<=(meanHeight+sigmaHeight)));
    Height_1sigma_Percent = length(Height_1sigma_Num) / length(endHeight);
    
    [mean_zExposed_fit, zExposed_stats] = robustfit(ones(length(zExposedAll),1),zExposedAll);
    mean_zExposed = max(0,mean_zExposed_fit(1));
    sigma_zExposed = zExposed_stats.s; % robust fit sigma
    mean_zExposed_ls = mean(zExposedAll); % ordinary mean
    sigma_zExposed_ls = zExposed_stats.ols_s; % rmse from least squares fit
        
    [mean_zDark_fit, zDark_stats] = robustfit(ones(length(zDarkAll),1),zDarkAll);
    mean_zDark = max(0,mean_zDark_fit(1));
    sigma_zDark = zDark_stats.s; % robust fit sigma
    mean_zDark_ls = mean(zDarkAll); % ordinary mean
    sigma_zDark_ls = zDark_stats.ols_s; % rmse from least squares fit
    
    % Average phase using "robustfit" to remove outlier
%     meanPhase2Pi_fit = robustfit(1:1:length(Phase2PiAll),Phase2PiAll);
    [meanPhase2Pi_fit, Phase2Pi_stats] = robustfit(ones(length(Phase2PiAll),1),Phase2PiAll);
    meanPhase2Pi = max(0,meanPhase2Pi_fit(1));
    sigmaPhase2Pi = Phase2Pi_stats.s; % robust fit sigma
    meanPhase2Pi_ls = mean(Phase2PiAll); % ordinary mean
    sigmaPhase2Pi_ls = Phase2Pi_stats.ols_s; % rmse from least squares fit
    
    % save results
    save(strcat(SaveFolder,'\MeasStats.mat'),'endHeight','meanHeight','meanHeight_SE','sigmaHeight','Height_3sigma_Num','Height_3sigma_Percent','Height_2sigma_Num','Height_2sigma_Percent','Height_1sigma_Num','Height_1sigma_Percent',...
        'meanHeight_ls','sigmaHeight_ls','meanHeight_fit','Height_stats',...
    'mean_zExposed','sigma_zExposed','mean_zExposed_ls','sigma_zExposed_ls','mean_zExposed_fit', 'zExposed_stats',...
    'mean_zDark','sigma_zDark','mean_zDark_ls','sigma_zDark_ls','mean_zDark_fit', 'zDark_stats',...
    'meanPhase2Pi','sigmaPhase2Pi','meanPhase2Pi_ls','sigmaPhase2Pi_ls','meanPhase2Pi_fit', 'Phase2Pi_stats');
end


if  isvector(m)
    figure(1) % Total height
    mNon0 = full(m); % nonzero elements in m
    mNon0(m==0) = NaN;
    plot(mNon0,'.');
    hold on
    plot(meanHeight*ones(length(mNon0),1),'r--');
    xlabel('Pixel Index'); ylabel('Cured Height (\mum)'); 
    legend('Cured Height of Each Pixel', sprintf('Average Height: %0.3f & Sigma: %0.3f', meanHeight,sigmaHeight));
    axis([0 length(mNon0)+5 0 round(max(HeightsAll(end,:)))+25]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_mean.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_mean.fig', SaveFolder));
    
    figure(2) % Exposed Height
    zExposedNon0 = full(zExposed); % nonzero elements in m
    zExposedNon0(zExposed==0) = NaN;
    plot(zExposedNon0,'.');
    hold on
    plot(mean_zExposed*ones(length(zExposedNon0),1),'r--');
    xlabel('Pixel Index'); ylabel('Cured Height under Exposure (\mum)'); 
    legend('Exposed Cured Height of Each Pixel', sprintf('Average Height: %0.3f & Sigma: %0.3f', mean_zExposed, sigma_zExposed));
    axis([0 length(zExposedNon0)+5 0 round(max(zExposed)+25)]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_mean.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_mean.fig', SaveFolder));
    
    figure(3) % Dark Height
    zDarkNon0 = full(zDark); % nonzero elements in m
    zDarkNon0(zDark==0) = NaN;
    plot(zDarkNon0,'.');
    hold on
    plot(mean_zDark*ones(length(zDarkNon0),1),'r--');
    xlabel('Pixel Index'); ylabel('Cured Height after Exposure (\mum)'); 
    legend('Dark Cured Height of Each Pixel', sprintf('Average Height: %0.3f & Sigma: %0.3f', mean_zDark, sigma_zDark));
    axis([0 length(zDarkNon0)+5 0 round(max(zDark)+25)]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_mean.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_mean.fig', SaveFolder));
    
    figure(4)% Total Height vs. Exposed Height
    plot(mNon0,'r.');
    hold on
    plot(meanHeight*ones(length(mNon0),1),'r--');
    hold on
    plot(zExposedNon0,'b*');
    hold on
    plot(mean_zExposed*ones(length(zExposedNon0),1),'b--');
    xlabel('Pixel Index'); ylabel('Cured Height (\mum)'); 
    legend('Total Cured Height of Each Pixel', sprintf('Total Cured Height with average of %0.3f', meanHeight),...
        'Exposed Cured Height of Each Pixel', sprintf('Exposed Cured Height with average of %0.3f', mean_zExposed));
    axis([0 length(mNon0)+5 0 round(max(HeightsAll(end,:)))+40]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed.fig', SaveFolder));
    
    figure(5) % Phase Angle
    phase2Pi_Non0 = full(phase2Pi); % nonzero elements in m
    phase2Pi_Non0(phase2Pi==0) = NaN;
    plot(phase2Pi_Non0,'X');
    hold on
    plot(meanPhase2Pi*ones(length(phase2Pi_Non0),1),'r--');
    xlabel('Pixel Index'); ylabel('Estimated Total Phase (unit:cycle, i.e. 2Pi rad)'); 
    legend('Pixel Phase Cycles', sprintf('Average Phase Cycles: %0.4f  & Sigma: %0.4f', meanPhase2Pi, sigmaPhase2Pi));
    axis([0 length(phase2Pi_Non0)+5 0 round(max(phase2Pi_Non0))+2]);
    saveas(gcf,sprintf('%s\\AllPointEstPhase_mean.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstPhase_mean.fig', SaveFolder));
else
    figure(1)% total height 3D view
    surf(WR, HR, m, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Cured Height (\mum)');
    zlim([0 round(max(HeightsAll(end,:)))+25]);
    colormap jet;
    colorbar;
    title({'Cured Height of Each Pixel', sprintf('Average Height: %0.3f (um) & Sigma: %0.3f (um)', meanHeight, sigmaHeight)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_3D.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_3D.fig', SaveFolder));

    figure(2) % total height 2D view
    surf(WR, HR, m, 'EdgeColor','None'), view(2), axis equal tight, colorbar;
    ylabel('Pixel Height Coordinate in Interferogram');xlabel('Pixel Width Coordinate in Interferogram'); zlabel('Cured Height (\mum)');
    colormap jet;
    title({'Cured Height of Each Pixel', sprintf('Average Height: %0.3f (um) & Sigma: %0.3f (um)', meanHeight, sigmaHeight)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_2D.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_2D.fig', SaveFolder));
    
    figure(3) % exposed cured height 3D view
    surf(WR, HR, zExposed, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Exposed Cured Height (\mum)');
    zlim([0 round(max(zExposedAll))+25]);
    colormap jet;
    colorbar;
    title({'Pixel Cured Height under Exposure', sprintf('Average Exposed Cured Height: %0.3f (um) & Sigma: %0.3f (um)', mean_zExposed, sigma_zExposed)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_3D.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_3D.fig', SaveFolder));
    
    figure(4) % dark cured height 3D view
    surf(WR, HR, zDark, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Dark Cured Height (\mum)');
    zlim([0 round(max(zDarkAll))+25]);
    colormap jet;
    colorbar;
    title({'Pixel Cured Height after Exposure', sprintf('Average Dark Cured Height: %0.3f (um) & Sigma: %0.3f (um)', mean_zDark, sigma_zDark)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_3D.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_3D.fig', SaveFolder));
    
    figure(5)% Total Height vs. Exposed Height in 3D view
    surf(WR, HR, m, 'EdgeColor','None');
    hold on,
    surf(WR, HR, zExposed, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Cured Height (\mum)');
    zlim([0 round(max(HeightsAll(end,:)))+25]);
    colormap jet;
    colorbar;
    title({'Pixel Cured Height: Total vs. Exposed',...
        sprintf('Top Graph: Total Height with average of %0.3f (um)', meanHeight),...
        sprintf('Bottom Graph: Exposed Cured Height with Average of %0.3f (um)', mean_zExposed)});
%     legend(sprintf('Total Height with average of %0.3f (um)', meanHeight),sprintf('Exposed Cured Height with Average of %0.3f (um)', mean_zExposed));  
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_3D.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_3D.fig', SaveFolder));
    
    figure(6) % phase angle 2D view
    surf(WR, HR, phase2Pi, 'EdgeColor','None'), view(2), axis equal tight, colorbar;
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Estimated Total Phase (unit:cycle, i.e. 2Pi rad)');
    colormap jet;
    title({'Pixel Phase Cycles', sprintf('Average Phase Cycles: %0.4f (unit:cycle, i.e. 2Pi rad)  & Sigma: %0.4f', meanPhase2Pi,sigmaPhase2Pi)});
    saveas(gcf,sprintf('%s\\AllPointEstPhase_2D.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstPhase_2D.fig', SaveFolder));
    
    

end


