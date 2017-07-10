function reportMeasureRet_Offline(MeasureRet,RunNo_uvClose)
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
SaveFolder = fullfile(ResultFolder,strcat('OfflineMeasRetPlots_',datestr(now,'yyyymmdd_HHMMSS')));
if ~exist(SaveFolder, 'dir')
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

    title({'Offline ICM&M Time Series of Grayscale',...
       sprintf('Pixel (Height,Width):(%03d, %03d)', PixelH, PixelW),...
       sprintf('Estimated Total Phase Angle: %0.4f Cycles; Cured Height: %0.3f um', Phase2PiAll(iPoint),endHeight(iPoint))});
    xlabel('Time(s)');ylabel('Grayscale');
    axis([0 max(dataX)+2 0 300]);

    pointYfig = sprintf('%s\\Grayscale_Px_%03d_%03d.png', SaveFolder, PixelH, PixelW);
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
%     plot(t, z,'Linewidth',2);
%     hold on
%     plot(t,z_uvClose*ones(length(t),1),'r--');
%     hold on
%     plot(t_uvClose*ones(length(z),1),z,'g--');
%     title({'Offline ICM&M Estimated Cured Height',...
%        sprintf('Pixel (Height,Width):(%03d, %03d)', PixelH, PixelW) });
%     xlabel('time (s)');ylabel('Cured Height (\mum)');
%     axis([0 max(TimesAll(end,:))+2 0 round(max(HeightsAll(end,:)))+20]);
%     legend(sprintf('Estimated Final Total Cured Height: %0.3f', z(end)),...
%        sprintf('Estimated Exposed Cured Height: %0.3f', z_uvClose) ,...
%        sprintf('Approximate Time when UV Lamp Closed: %0.3f', t_uvClose));
%     pointZfig = sprintf('%s\\ICMM_Height_Px_%03d_%03d.png', SaveFolder, PixelH, PixelW);
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
%     plot(t, [0;freq],'Linewidth',2);
%     title({'Offline ICM&M Estimated Instantaneous Frequency (Hz)',...
%        sprintf('Pixel (Height,Width):(%03d, %03d)', PixelH, PixelW) });
%     xlabel('time (s)');ylabel('Instantaneous Frequency (Hz)');
%     axis([0 max(TimesAll(end,:))+2 0 2]);
%     legend(sprintf('Estimated Total Phase Angle: %0.4f Cycles', Phase2PiAll(iPoint)));
%     pointFfig = sprintf('%s\\ICMM_Freq_Px_%03d_%03d.png', SaveFolder, PixelH, PixelW);
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
filteredHeightProfile = zeros(lenH, lenW);
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

%% 2-D median filtering across the ROI
filteredHeightProfile = medfilt2(m);
filtered_zExposed = medfilt2(zExposed);
filtered_zDark = medfilt2(zDark);
filtered_phase2Pi = medfilt2(phase2Pi);
% %%
% m = filteredHeightProfile;
% zExposed = filtered_zExposed;
% zDark = filtered_zDark;
% phase2Pi = filtered_phase2Pi;

%% Calculate and save the statistics of initial ICM&M results: to estimate deviation for outlier detection
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
    
    save(strcat(SaveFolder,'\InitialMeasStats.mat'), 'endHeight','meanHeight_ls','sigmaHeight_ls',...
    'mean_zExposed_ls','sigma_zExposed_ls',...
    'mean_zDark_ls','sigma_zDark_ls',...
    'meanPhase2Pi_ls','sigmaPhase2Pi_ls');

    % for plot, to be consistent with robustfit results plot, use same terms
    meanHeight = meanHeight_ls;
    final_sigmaHeight = sigmaHeight_ls;    
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
    save(strcat(SaveFolder,'\InitialMeasStats.mat'),'endHeight','meanHeight','meanHeight_SE','sigmaHeight','Height_3sigma_Num','Height_3sigma_Percent','Height_2sigma_Num','Height_2sigma_Percent','Height_1sigma_Num','Height_1sigma_Percent',...
        'meanHeight_ls','sigmaHeight_ls','meanHeight_fit','Height_stats',...
    'mean_zExposed','sigma_zExposed','mean_zExposed_ls','sigma_zExposed_ls','mean_zExposed_fit', 'zExposed_stats',...
    'mean_zDark','sigma_zDark','mean_zDark_ls','sigma_zDark_ls','mean_zDark_fit', 'zDark_stats',...
    'meanPhase2Pi','sigmaPhase2Pi','meanPhase2Pi_ls','sigmaPhase2Pi_ls','meanPhase2Pi_fit', 'Phase2Pi_stats');
end

%% Outlier Detection and Treatment
finalHeightProfile = m;
final_zExposed = zExposed;
final_zDark = zDark;
final_phase2Pi = phase2Pi;
for jH = 1:1:lenH
    for jW = 1:1:lenW
        if (m(jH,jW)> meanHeight + sigmaHeight) || (m(jH,jW)< meanHeight - sigmaHeight)  % Outlier detection: off 1-sigma
            finalHeightProfile(jH,jW) = filteredHeightProfile(jH,jW); % Outlier treatment: replace with filtered result
            final_zExposed(jH,jW) = filtered_zExposed(jH,jW);
            final_zDark(jH,jW) = filtered_zDark(jH,jW);
            final_phase2Pi(jH,jW) = filtered_phase2Pi(jH,jW);
        end
    end
end
InitialHeightProfile = m;
Initial_zExposed = zExposed;
Initial_zDark = zDark;
Initial_phase2Pi = phase2Pi;
    save(strcat(SaveFolder,'\ProfileMeasResults.mat'), 'InitialHeightProfile','finalHeightProfile','filteredHeightProfile',...
   'Initial_zExposed','final_zExposed','filtered_zExposed',...
    'Initial_zDark','final_zDark','filtered_zDark',...
    'Initial_phase2Pi','final_phase2Pi','filtered_phase2Pi');

%% Statistics of the final results
if length(MeasureRet)<3
    % Average height
%     endHeight = HeightsAll(end,:);
    final_meanHeight_ls = mean(finalHeightProfile(:));
    final_sigmaHeight_ls = std(finalHeightProfile(:)); % rmse from least squares fit
    
    final_mean_zExposed_ls = mean(final_zExposed(:));
    final_sigma_zExposed_ls = std(final_zExposed(:)); % rmse from least squares fit
    
    final_mean_zDark_ls = mean(final_zDark(:));
    final_sigma_zDark_ls = std(final_zDark(:)); % rmse from least squares fit
    
    % Average phase
    final_meanPhase2Pi_ls = mean(final_phase2Pi(:));
    final_sigmaPhase2Pi_ls = std(final_phase2Pi(:)); % rmse from least squares fit
    
    save(strcat(SaveFolder,'\FinalMeasStats.mat'), 'finalHeightProfile','final_meanHeight_ls','final_sigmaHeight_ls',...
    'final_mean_zExposed_ls','final_sigma_zExposed_ls',...
    'final_mean_zDark_ls','final_sigma_zDark_ls',...
    'final_meanPhase2Pi_ls','final_sigmaPhase2Pi_ls');

    % for plot, to be consistent with robustfit results plot, use same terms
    final_meanHeight = final_meanHeight_ls;
    final_sigmaHeight = final_sigmaHeight_ls;    
    final_mean_zExposed = final_mean_zExposed_ls;
    final_sigma_zExposed = final_sigma_zExposed_ls;     
    final_mean_zDark = final_mean_zDark_ls;
    final_sigma_zDark = final_sigma_zDark_ls;     
    final_meanPhase2Pi = final_meanPhase2Pi_ls;
    final_sigmaPhase2Pi = final_sigmaPhase2Pi_ls;
else
    % Average height using "robustfit" to remove outlier
%     meanHeight_fit = robustfit(1:1:length(HeightsAll(end,:)),HeightsAll(end,:));
    [final_meanHeight_fit,final_Height_stats] = robustfit(ones(length(finalHeightProfile(:)),1),finalHeightProfile(:));
    final_meanHeight = max(0,final_meanHeight_fit(1)); % robust fit mean
    final_meanHeight_SE = final_Height_stats.se(1); % standard error of meanHeight estimates
    final_sigmaHeight = final_Height_stats.s; % robust fit sigma
    final_meanHeight_ls = mean(finalHeightProfile(:)); % ordinary mean
    final_sigmaHeight_ls = final_Height_stats.ols_s; % rmse from least squares fit
    endHeight = finalHeightProfile(:);
    final_Height_3sigma_Num=endHeight(find(endHeight>=(final_meanHeight-3*final_sigmaHeight)& endHeight<=(final_meanHeight+3*final_sigmaHeight)));
    final_Height_3sigma_Percent = length(final_Height_3sigma_Num) / length(endHeight);
    final_Height_2sigma_Num=endHeight(find(endHeight>=(final_meanHeight-2*final_sigmaHeight)& endHeight<=(final_meanHeight+2*final_sigmaHeight)));
    final_Height_2sigma_Percent = length(final_Height_2sigma_Num) / length(endHeight);
    final_Height_1sigma_Num=endHeight(find(endHeight>=(final_meanHeight-final_sigmaHeight)& endHeight<=(final_meanHeight+final_sigmaHeight)));
    final_Height_1sigma_Percent = length(final_Height_1sigma_Num) / length(endHeight);
    
    [final_mean_zExposed_fit, final_zExposed_stats] = robustfit(ones(length(final_zExposed(:)),1),final_zExposed(:));
    final_mean_zExposed = max(0,final_mean_zExposed_fit(1));
    final_sigma_zExposed = final_zExposed_stats.s; % robust fit sigma
    final_mean_zExposed_ls = mean(final_zExposed(:)); % ordinary mean
    final_sigma_zExposed_ls = final_zExposed_stats.ols_s; % rmse from least squares fit
        
    [final_mean_zDark_fit, final_zDark_stats] = robustfit(ones(length(final_zDark(:)),1),final_zDark(:));
    final_mean_zDark = max(0,final_mean_zDark_fit(1));
    final_sigma_zDark = final_zDark_stats.s; % robust fit sigma
    final_mean_zDark_ls = mean(final_zDark(:)); % ordinary mean
    final_sigma_zDark_ls = final_zDark_stats.ols_s; % rmse from least squares fit
    
    % Average phase using "robustfit" to remove outlier
%     meanPhase2Pi_fit = robustfit(1:1:length(Phase2PiAll),Phase2PiAll);
    [final_meanPhase2Pi_fit, final_Phase2Pi_stats] = robustfit(ones(length(final_phase2Pi(:)),1),final_phase2Pi(:));
    final_meanPhase2Pi = max(0,final_meanPhase2Pi_fit(1));
    final_sigmaPhase2Pi = final_Phase2Pi_stats.s; % robust fit sigma
    final_meanPhase2Pi_ls = mean(final_phase2Pi(:)); % ordinary mean
    final_sigmaPhase2Pi_ls = final_Phase2Pi_stats.ols_s; % rmse from least squares fit
    
    % save results
    save(strcat(SaveFolder,'\FinalMeasStats.mat'),'finalHeightProfile','final_meanHeight','final_meanHeight_SE','final_sigmaHeight','final_Height_3sigma_Num','final_Height_3sigma_Percent','final_Height_2sigma_Num','final_Height_2sigma_Percent','final_Height_1sigma_Num','final_Height_1sigma_Percent',...
        'final_meanHeight_ls','final_sigmaHeight_ls','final_meanHeight_fit','final_Height_stats',...
    'final_mean_zExposed','final_sigma_zExposed','final_mean_zExposed_ls','final_sigma_zExposed_ls','final_mean_zExposed_fit', 'final_zExposed_stats',...
    'final_mean_zDark','final_sigma_zDark','final_mean_zDark_ls','final_sigma_zDark_ls','final_mean_zDark_fit', 'final_zDark_stats',...
    'final_meanPhase2Pi','final_sigmaPhase2Pi','final_meanPhase2Pi_ls','final_sigmaPhase2Pi_ls','final_meanPhase2Pi_fit', 'final_Phase2Pi_stats');
end

%% Plots initial ICM&M results
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
    saveas(gcf,sprintf('%s\\AllPointEstHeights_mean_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_mean_Initial.fig', SaveFolder));
    
    figure(2) % Exposed Height
    zExposedNon0 = full(zExposed); % nonzero elements in m
    zExposedNon0(zExposed==0) = NaN;
    plot(zExposedNon0,'.');
    hold on
    plot(mean_zExposed*ones(length(zExposedNon0),1),'r--');
    xlabel('Pixel Index'); ylabel('Cured Height under Exposure (\mum)'); 
    legend('Exposed Cured Height of Each Pixel', sprintf('Average Height: %0.3f & Sigma: %0.3f', mean_zExposed, sigma_zExposed));
    axis([0 length(zExposedNon0)+5 0 round(max(zExposed)+25)]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_mean_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_mean_Initial.fig', SaveFolder));
    
    figure(3) % Dark Height
    zDarkNon0 = full(zDark); % nonzero elements in m
    zDarkNon0(zDark==0) = NaN;
    plot(zDarkNon0,'.');
    hold on
    plot(mean_zDark*ones(length(zDarkNon0),1),'r--');
    xlabel('Pixel Index'); ylabel('Cured Height after Exposure (\mum)'); 
    legend('Dark Cured Height of Each Pixel', sprintf('Average Height: %0.3f & Sigma: %0.3f', mean_zDark, sigma_zDark));
    axis([0 length(zDarkNon0)+5 0 round(max(zDark)+25)]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_mean_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_mean_Initial.fig', SaveFolder));
    
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
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_Initial.fig', SaveFolder));
    
    figure(5) % Phase Angle
    phase2Pi_Non0 = full(phase2Pi); % nonzero elements in m
    phase2Pi_Non0(phase2Pi==0) = NaN;
    plot(phase2Pi_Non0,'X');
    hold on
    plot(meanPhase2Pi*ones(length(phase2Pi_Non0),1),'r--');
    xlabel('Pixel Index'); ylabel('Estimated Total Phase (unit:cycle, i.e. 2Pi rad)'); 
    legend('Pixel Phase Cycles', sprintf('Average Phase Cycles: %0.4f  & Sigma: %0.4f', meanPhase2Pi, sigmaPhase2Pi));
    axis([0 length(phase2Pi_Non0)+5 0 round(max(phase2Pi_Non0))+2]);
    saveas(gcf,sprintf('%s\\AllPointEstPhase_mean_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstPhase_mean_Initial.fig', SaveFolder));
else
    figure(1)% total height 3D view
    surf(WR, HR, m, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Cured Height (\mum)');
    zlim([0 round(max(HeightsAll(end,:)))+25]);
    colormap jet;
    colorbar;
    title({'Cured Height of Each Pixel', sprintf('Average Height: %0.3f (um) & Sigma: %0.3f (um)', meanHeight, sigmaHeight)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_3D_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_3D_Initial.fig', SaveFolder));

    figure(2) % total height 2D view
    surf(WR, HR, m, 'EdgeColor','None'), view(2), axis equal tight, colorbar;
    ylabel('Pixel Height Coordinate in Interferogram');xlabel('Pixel Width Coordinate in Interferogram'); zlabel('Cured Height (\mum)');
    colormap jet;
    title({'Cured Height of Each Pixel', sprintf('Average Height: %0.3f (um) & Sigma: %0.3f (um)', meanHeight, sigmaHeight)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_2D_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_2D_Initial.fig', SaveFolder));
    
    figure(3) % exposed cured height 3D view
    surf(WR, HR, zExposed, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Exposed Cured Height (\mum)');
    zlim([0 round(max(zExposedAll))+25]);
    colormap jet;
    colorbar;
    title({'Pixel Cured Height under Exposure', sprintf('Average Exposed Cured Height: %0.3f (um) & Sigma: %0.3f (um)', mean_zExposed, sigma_zExposed)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_3D_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_3D_Initial.fig', SaveFolder));
    
    figure(4) % dark cured height 3D view
    surf(WR, HR, zDark, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Dark Cured Height (\mum)');
    zlim([0 round(max(zDarkAll))+25]);
    colormap jet;
    colorbar;
    title({'Pixel Cured Height after Exposure', sprintf('Average Dark Cured Height: %0.3f (um) & Sigma: %0.3f (um)', mean_zDark, sigma_zDark)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_3D_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_3D_Initial.fig', SaveFolder));
    
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
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_3D_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_3D_Initial.fig', SaveFolder));
    
    figure(6) % phase angle 2D view
    surf(WR, HR, phase2Pi, 'EdgeColor','None'), view(2), axis equal tight, colorbar;
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Estimated Total Phase (unit:cycle, i.e. 2Pi rad)');
    colormap jet;
    title({'Pixel Phase Cycles', sprintf('Average Phase Cycles: %0.4f (unit:cycle, i.e. 2Pi rad)  & Sigma: %0.4f', meanPhase2Pi,sigmaPhase2Pi)});
    saveas(gcf,sprintf('%s\\AllPointEstPhase_2D_Initial.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstPhase_2D_Initial.fig', SaveFolder));
end
%% Plots final ICM&M results
if  isvector(finalHeightProfile)
    figure(6) % Total height
    final_mNon0 = full(finalHeightProfile); % nonzero elements in m
    final_mNon0(finalHeightProfile==0) = NaN;
    plot(final_mNon0,'.');
    hold on
    plot(final_meanHeight*ones(length(final_mNon0),1),'r--');
    xlabel('Pixel Index'); ylabel('Cured Height (\mum)'); 
    legend('Cured Height of Each Pixel', sprintf('Average Height: %0.3f & Sigma: %0.3f', final_meanHeight,final_sigmaHeight));
    axis([0 length(final_mNon0)+5 0 round(max(finalHeightProfile(:)))+25]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_mean_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_mean_Final.fig', SaveFolder));
    
    figure(7) % Exposed Height
    final_zExposedNon0 = full(final_zExposed); % nonzero elements in m
    final_zExposedNon0(final_zExposed==0) = NaN;
    plot(final_zExposedNon0,'.');
    hold on
    plot(final_mean_zExposed*ones(length(final_zExposedNon0),1),'r--');
    xlabel('Pixel Index'); ylabel('Cured Height under Exposure (\mum)'); 
    legend('Exposed Cured Height of Each Pixel', sprintf('Average Height: %0.3f & Sigma: %0.3f', final_mean_zExposed, final_sigma_zExposed));
    axis([0 length(final_zExposedNon0)+5 0 round(max(final_zExposed(:))+25)]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_mean_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_mean_Final.fig', SaveFolder));
    
    figure(8) % Dark Height
    final_zDarkNon0 = full(final_zDark); % nonzero elements in m
    final_zDarkNon0(final_zDark==0) = NaN;
    plot(final_zDarkNon0,'.');
    hold on
    plot(final_mean_zDark*ones(length(final_zDarkNon0),1),'r--');
    xlabel('Pixel Index'); ylabel('Cured Height after Exposure (\mum)'); 
    legend('Dark Cured Height of Each Pixel', sprintf('Average Height: %0.3f & Sigma: %0.3f', final_mean_zDark, final_sigma_zDark));
    axis([0 length(final_zDarkNon0)+5 0 round(max(final_zDark(:))+25)]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_mean_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_mean_Final.fig', SaveFolder));
    
    figure(9)% Total Height vs. Exposed Height
    plot(final_mNon0,'r.');
    hold on
    plot(final_meanHeight*ones(length(final_mNon0),1),'r--');
    hold on
    plot(final_zExposedNon0,'b*');
    hold on
    plot(final_mean_zExposed*ones(length(final_zExposedNon0),1),'b--');
    xlabel('Pixel Index'); ylabel('Cured Height (\mum)'); 
    legend('Total Cured Height of Each Pixel', sprintf('Total Cured Height with average of %0.3f', final_meanHeight),...
        'Exposed Cured Height of Each Pixel', sprintf('Exposed Cured Height with average of %0.3f', final_mean_zExposed));
    axis([0 length(final_mNon0)+5 0 round(max(finalHeightProfile(:)))+40]);
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_Final.fig', SaveFolder));
    
    figure(10) % Phase Angle
    final_phase2Pi_Non0 = full(final_phase2Pi); % nonzero elements in m
    final_phase2Pi_Non0(final_phase2Pi==0) = NaN;
    plot(final_phase2Pi_Non0,'X');
    hold on
    plot(final_meanPhase2Pi*ones(length(final_phase2Pi_Non0),1),'r--');
    xlabel('Pixel Index'); ylabel('Estimated Total Phase (unit:cycle, i.e. 2Pi rad)'); 
    legend('Pixel Phase Cycles', sprintf('Average Phase Cycles: %0.4f  & Sigma: %0.4f', final_meanPhase2Pi, final_sigmaPhase2Pi));
    axis([0 length(final_phase2Pi_Non0)+5 0 round(max(final_phase2Pi_Non0(:)))+2]);
    saveas(gcf,sprintf('%s\\AllPointEstPhase_mean_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstPhase_mean_Final.fig', SaveFolder));
else
    figure(7)% total height 3D view
    surf(WR, HR, finalHeightProfile, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Cured Height (\mum)');
    zlim([0 round(max(finalHeightProfile(:)))+25]);
    colormap jet;
    colorbar;
    title({'Cured Height of Each Pixel', sprintf('Average Height: %0.3f (um) & Sigma: %0.3f (um)', final_meanHeight, final_sigmaHeight)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_3D_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_3D_Final.fig', SaveFolder));

    figure(8) % total height 2D view
    surf(WR, HR, finalHeightProfile, 'EdgeColor','None'), view(2), axis equal tight, colorbar;
    ylabel('Pixel Height Coordinate in Interferogram');xlabel('Pixel Width Coordinate in Interferogram'); zlabel('Cured Height (\mum)');
    colormap jet;
    title({'Cured Height of Each Pixel', sprintf('Average Height: %0.3f (um) & Sigma: %0.3f (um)', final_meanHeight, final_sigmaHeight)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_2D_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_2D_Final.fig', SaveFolder));
    
    figure(9) % exposed cured height 3D view
    surf(WR, HR, final_zExposed, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Exposed Cured Height (\mum)');
    zlim([0 round(max(final_zExposed(:)))+25]);
    colormap jet;
    colorbar;
    title({'Pixel Cured Height under Exposure', sprintf('Average Exposed Cured Height: %0.3f (um) & Sigma: %0.3f (um)', final_mean_zExposed, final_sigma_zExposed)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_3D_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zExposed_3D_Final.fig', SaveFolder));
    
    figure(10) % dark cured height 3D view
    surf(WR, HR, final_zDark, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Dark Cured Height (\mum)');
    zlim([0 round(max(final_zDark(:)))+25]);
    colormap jet;
    colorbar;
    title({'Pixel Cured Height after Exposure', sprintf('Average Dark Cured Height: %0.3f (um) & Sigma: %0.3f (um)', final_mean_zDark, final_sigma_zDark)});
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_3D_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_zDark_3D_Final.fig', SaveFolder));
    
    figure(11)% Total Height vs. Exposed Height in 3D view
    surf(WR, HR, finalHeightProfile, 'EdgeColor','None');
    hold on,
    surf(WR, HR, final_zExposed, 'EdgeColor','None');
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Cured Height (\mum)');
    zlim([0 round(max(finalHeightProfile(:)))+25]);
    colormap jet;
    colorbar;
    title({'Pixel Cured Height: Total vs. Exposed',...
        sprintf('Top Graph: Total Height with average of %0.3f (um)', final_meanHeight),...
        sprintf('Bottom Graph: Exposed Cured Height with Average of %0.3f (um)', final_mean_zExposed)});
%     legend(sprintf('Total Height with average of %0.3f (um)', meanHeight),sprintf('Exposed Cured Height with Average of %0.3f (um)', mean_zExposed));  
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_3D_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstHeights_Total_Exposed_3D_Final.fig', SaveFolder));
    
    figure(12) % phase angle 2D view
    surf(WR, HR, final_phase2Pi, 'EdgeColor','None'), view(2), axis equal tight, colorbar;
    ylabel('Pixel Height Coord');xlabel('Pixel Width Coord'); zlabel('Estimated Total Phase (unit:cycle, i.e. 2Pi rad)');
    colormap jet;
    title({'Pixel Phase Cycles', sprintf('Average Phase Cycles: %0.4f (unit:cycle, i.e. 2Pi rad)  & Sigma: %0.4f', final_meanPhase2Pi,final_sigmaPhase2Pi)});
    saveas(gcf,sprintf('%s\\AllPointEstPhase_2D_Final.png', SaveFolder));
    saveas(gcf,sprintf('%s\\AllPointEstPhase_2D_Final.fig', SaveFolder));
end



