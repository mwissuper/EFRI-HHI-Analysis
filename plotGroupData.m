%% Plot results for all subj's combined

clear; close all; clc;
load HHI2017_LateStats_force_MW.mat;
% Find means for each trial type for each subj
n = 0;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)
conds = {'Solo Ground','Assist Ground','Solo Beam','Assist Beam'};

%% Concatenate force, vel, power data: do only for subj's with good force data and for assist cond's

n = 0;
for subj = subj_array_force
    n = n + 1;
    for i = [2 4]
        rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds{i});
        %% Interaction point X dir
        % Power
        powerAbsIntPtX(n,i/2) = nanmean(LateGroup.meanAbsPowerIntPtX(rows));
        SDpowerAbsIntPtX(n,i/2) = nanmean(LateGroup.SDAbsPowerIntPtX(rows));
        powerPosIntPtX(n,i/2) = nanmean(LateGroup.meanPosPowerIntPtX(rows));
        powerNegIntPtX(n,i/2) = nanmean(LateGroup.meanNegPowerIntPtX(rows));
%         % Sum of work for trial (biased by trial length)
%         PosWorkIntPtX(n,i/2) = nanmean(LateGroup.PosWorkIntPtX(rows));
%         NegWorkIntPtX(n,i/2) = nanmean(LateGroup.NegWorkIntPtX(rows));
%         % Mean of work vs. time for trial
%         meanWorkIntPtX(n,i/2) = nanmean(LateGroup.meanWorkIntPtX(rows)); % Net work
%         meanWorkMagIntPtX(n,i/2) = nanmean(LateGroup.meanWorkMagIntPtX(rows)); % Abs work
%         meanPosWorkIntPtX(n,i/2) = nanmean(LateGroup.meanPosWorkIntPtX(rows));
%         meanNegWorkIntPtX(n,i/2) = nanmean(LateGroup.meanNegWorkIntPtX(rows));
        % Mean force vs. time for trial
        meanFx(n,i/2) = nanmean(LateGroup.meanFx(rows));
        SDFx(n,i/2) = nanmean(LateGroup.SDFx(rows));
        meanPosFx(n,i/2) = nanmean(LateGroup.meanPosFx(rows));
        meanNegFx(n,i/2) = nanmean(LateGroup.meanNegFx(rows));
        % SD force vs. time for trial
        SDFx(n,i/2) = nanstd(LateGroup.SDFx(rows));
        SDPosFx(n,i/2) = nanstd(LateGroup.SDPosFx(rows));
        SDNegFx(n,i/2) = nanstd(LateGroup.SDNegFx(rows));
        % Velocity R FIN marker
        meanVx(n,i/2) = nanmean(LateGroup.meanVx(rows));
        meanPosVx(n,i/2) = nanmean(LateGroup.meanPosVx(rows));
        meanNegVx(n,i/2) = nanmean(LateGroup.meanNegVx(rows));
        % corr power to Fresid
        meanrPowerFresX(n,i/2) = nanmean(LateGroup.rPowerFresX(rows));
        % xcorr IP to POB Clav
        meanxcorrFIPClavX(n,i/2) = nanmean(LateGroup.xcorrFIPClavX(rows));
        meanxcorrFIPvClavX(n,i/2) = nanmean(LateGroup.xcorrFIPvClavX(rows));
        meanxcorrvIPvClavX(n,i/2) = nanmean(LateGroup.xcorrvIPvClavX(rows));
        % xcorr lag IP to POB Clav
        meanLagFIPClavX(n,i/2) = nanmean(LateGroup.lagFIPClavX(rows));
        meanLagFIPvClavX(n,i/2) = nanmean(LateGroup.lagFIPvClavX(rows));
        meanLagvIPvClavX(n,i/2) = nanmean(LateGroup.lagvIPvClavX(rows));
        % POB armLen (not force but related to xcorr above)
        meanArmPOBX(n,i/2) = nanmean(LateGroup.meanArmPOBX(rows));
        SDarmPOBX(n,i/2) = nanmean(LateGroup.SDarmPOBX(rows));
        
        %% Interaction point Y dir
        % Power
        powerAbsIntPtY(n,i/2) = nanmean(LateGroup.meanAbsPowerIntPtY(rows));
        powerPosIntPtY(n,i/2) = nanmean(LateGroup.meanPosPowerIntPtY(rows));
        powerNegIntPtY(n,i/2) = nanmean(LateGroup.meanNegPowerIntPtY(rows));
%         % Sum of work for trial (biased by trial length)
%         PosWorkIntPtY(n,i/2) = nanmean(LateGroup.PosWorkIntPtY(rows));
%         NegWorkIntPtY(n,i/2) = nanmean(LateGroup.NegWorkIntPtY(rows));
%         % Mean of work vs. time for trial
%         meanWorkIntPtY(n,i/2) = nanmean(LateGroup.meanWorkIntPtY(rows)); % Net work
%         meanWorkMagIntPtY(n,i/2) = nanmean(LateGroup.meanWorkMagIntPtY(rows)); % Abs work
%         meanPosWorkIntPtY(n,i/2) = nanmean(LateGroup.meanPosWorkIntPtY(rows));
%         meanNegWorkIntPtY(n,i/2) = nanmean(LateGroup.meanNegWorkIntPtY(rows));
        % Mean force vs. time for trial
        meanFy(n,i/2) = nanmean(LateGroup.meanFy(rows));
        meanPosFy(n,i/2) = nanmean(LateGroup.meanPosFy(rows));
        meanNegFy(n,i/2) = nanmean(LateGroup.meanNegFy(rows));
        % SD force vs. time for trial
        SDFy(n,i/2) = nanstd(LateGroup.SDFy(rows));
        SDPosFy(n,i/2) = nanstd(LateGroup.SDPosFy(rows));
        SDNegFy(n,i/2) = nanstd(LateGroup.SDNegFy(rows));
        % Velocity R FIN marker
        meanVy(n,i/2) = nanmean(LateGroup.meanVy(rows));
        meanPosVy(n,i/2) = nanmean(LateGroup.meanPosVy(rows));
        meanNegVy(n,i/2) = nanmean(LateGroup.meanNegVy(rows));
        
        %% Interaction point Z dir
        % Power
        powerAbsIntPtZ(n,i/2) = nanmean(LateGroup.meanAbsPowerIntPtZ(rows));
        powerPosIntPtZ(n,i/2) = nanmean(LateGroup.meanPosPowerIntPtZ(rows));
        powerNegIntPtZ(n,i/2) = nanmean(LateGroup.meanNegPowerIntPtZ(rows));
%         % Sum of work for trial (biased by trial length)
%         PosWorkIntPtZ(n,i/2) = nanmean(LateGroup.PosWorkIntPtZ(rows));
%         NegWorkIntPtZ(n,i/2) = nanmean(LateGroup.NegWorkIntPtZ(rows));
%         % Mean of work vs. time for trial
%         meanWorkIntPtZ(n,i/2) = nanmean(LateGroup.meanWorkIntPtZ(rows)); % Net work
%         meanWorkMagIntPtZ(n,i/2) = nanmean(LateGroup.meanWorkMagIntPtZ(rows)); % Abs work
%         meanPosWorkIntPtZ(n,i/2) = nanmean(LateGroup.meanPosWorkIntPtZ(rows));
%         meanNegWorkIntPtZ(n,i/2) = nanmean(LateGroup.meanNegWorkIntPtZ(rows));
        % Mean force vs. time for trial
        meanFz(n,i/2) = nanmean(LateGroup.meanFz(rows));
        meanPosFz(n,i/2) = nanmean(LateGroup.meanPosFz(rows));
        meanNegFz(n,i/2) = nanmean(LateGroup.meanNegFz(rows));
        % SD force vs. time for trial
        SDFz(n,i/2) = nanstd(LateGroup.SDFz(rows));
        SDPosFz(n,i/2) = nanstd(LateGroup.SDPosFz(rows));
        SDNegFz(n,i/2) = nanstd(LateGroup.SDNegFz(rows));
        % Velocity R FIN marker
        meanVz(n,i/2) = nanmean(LateGroup.meanVz(rows));
        meanPosVz(n,i/2) = nanmean(LateGroup.meanPosVz(rows));
        meanNegVz(n,i/2) = nanmean(LateGroup.meanNegVz(rows));

        %% Clavicle x/ML dir
        % Velocity CLAV marker
        meanVxPOB(n,i/2) = nanmean(LateGroup.meanVx_clav(rows));
        % Regression rFin to force (nan value if n.s. fit)
        meanMxPOB(n,i/2) = nanmean(LateGroup.mx_clav(rows));
        meanBxPOB(n,i/2) = nanmean(LateGroup.bx_clav(rows));
        meanKxPOB(n,i/2) = nanmean(LateGroup.kx_clav(rows));
        meanRsqxPOB(n,i/2) = nanmean(LateGroup.Rsqx_clav(rows));
        
        %% Clavicle y/AP dir
        % Velocity CLAV marker
        meanVyPOB(n,i/2) = nanmean(LateGroup.meanVy_clav(rows));
        % Regression rFin to force (nan value if n.s. fit)
        meanMyPOB(n,i/2) = nanmean(LateGroup.my_clav(rows));
        meanByPOB(n,i/2) = nanmean(LateGroup.by_clav(rows));
        meanKyPOB(n,i/2) = nanmean(LateGroup.ky_clav(rows));
        meanRsqyPOB(n,i/2) = nanmean(LateGroup.Rsqy_clav(rows));
        
        %% Clavicle z/Vert. dir
        % Velocity CLAV marker
        meanVzPOB(n,i/2) = nanmean(LateGroup.meanVz_clav(rows));
        % Regression rFin to force (nan value if n.s. fit)
        meanMzPOB(n,i/2) = nanmean(LateGroup.mz_clav(rows));
        meanBzPOB(n,i/2) = nanmean(LateGroup.bz_clav(rows));
        meanKzPOB(n,i/2) = nanmean(LateGroup.kz_clav(rows));
        meanRsqzPOB(n,i/2) = nanmean(LateGroup.Rsqz_clav(rows));
        
        %% Calculate percent of trials where param was sig.
        if i == 4 % Assist Beam condition only
            mxPOBsig(n) = length(find(~isnan(LateGroup.mx_clav(rows))))/length(LateGroup.mx_clav(rows));
            bxPOBsig(n) = length(find(~isnan(LateGroup.bx_clav(rows))))/length(LateGroup.bx_clav(rows));
            kxPOBsig(n) = length(find(~isnan(LateGroup.kx_clav(rows))))/length(LateGroup.kx_clav(rows));
            mzPOBsig(n) = length(find(~isnan(LateGroup.mz_clav(rows))))/length(LateGroup.mz_clav(rows));
            bzPOBsig(n) = length(find(~isnan(LateGroup.bz_clav(rows))))/length(LateGroup.bz_clav(rows));
            kzPOBsig(n) = length(find(~isnan(LateGroup.kz_clav(rows))))/length(LateGroup.kz_clav(rows));
        end
    end
end

save LateStats_groupMeans_force

%% Plot regression POB R^2 and coeff values ea dir Partner Beam
plotind = 0;

% ML dir
plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanRsqxPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Rsq'); 
title('Fit POB clav ML direction');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanMxPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Mass (kg)'); hline(0,'k');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanBxPOB(:,2)); set(gca,'xticklabel','All partnerships');box off; set(gca,'tickdir','out');
ylabel('Damping (Ns/m)'); hline(0,'k');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanKxPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Stiffness (N/m)'); hline(0,'k');

% AP dir
plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanRsqyPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Rsq'); 
title('Fit POB clav AP direction');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanMyPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Mass (kg)'); hline(0,'k');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanByPOB(:,2)); set(gca,'xticklabel','All partnerships');box off; set(gca,'tickdir','out');
ylabel('Damping (Ns/m)'); hline(0,'k');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanKyPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Stiffness (N/m)'); hline(0,'k');

% Vert dir
plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanRsqzPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Rsq');
title('Fit POB clav vert direction');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanMzPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Mass (kg)'); hline(0,'k');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanBzPOB(:,2)); set(gca,'xticklabel','All partnerships');box off; set(gca,'tickdir','out');
ylabel('Damping (Ns/m)'); hline(0,'k');

plotind = plotind + 1;
subplot(3,4,plotind)
boxplot(meanKzPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('Stiffness (N/m)'); hline(0,'k');

%% mean percent Partner Beam trials where m, b, k sig. ML dir

plotind = 0;
numrows = 3; numcols = 1;

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
bar(mxPOBsig*100),ylabel('% trials m sig')
box off; ylim([0 100]); xlabel('Partnership');
set(gca,'tickdir','out');
title('Fit POB clav');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
bar(bxPOBsig*100),ylabel('% trials b sig')
box off; ylim([0 100]); xlabel('Partnership');
set(gca,'tickdir','out');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
bar(kxPOBsig*100),ylabel('% trials k sig')
box off; ylim([0 100]); xlabel('Partnership');
set(gca,'tickdir','out');

%% Plot correlation power to regression residual for ML dir

boxplot(meanrPowerFresX(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
ylabel('r');
title('Correlation between ML power and F model residual');

%% Plot xcorr lag different IP to POB Clav signals for ML dir
xtlab{1} = conds{2};
xtlab{2} = conds{4};

subplot(2,3,1)
boxplot(meanLagFIPClavX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');
ylabel('lag at max xcorr (s)');
title('ML force and Clav pos');

subplot(2,3,2)
boxplot(meanLagFIPvClavX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');
title('ML force and Clav vel');

subplot(2,3,3)
boxplot(meanLagvIPvClavX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');
title('IP and Clav vels');

subplot(2,3,4)
boxplot(meanxcorrFIPClavX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');
ylabel('max xcorr');

subplot(2,3,5)
boxplot(meanxcorrFIPvClavX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');

subplot(2,3,6)
boxplot(meanxcorrvIPvClavX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');

%% Plot POB arm len to check xcorr results
xtlab{1} = conds{2};
xtlab{2} = conds{4};

subplot(1,2,1)
boxplot(meanArmPOBX),set(gca,'xticklabel',xtlab),hold on
box off; set(gca,'tickdir','out'),title('Distance IP to POB Clav (m)');
ylabel('mean (m)');

subplot(1,2,2)
boxplot(SDarmPOBX),set(gca,'xticklabel',xtlab),hold on
box off; set(gca,'tickdir','out')
ylabel('SD (m)');

%% Plot mean abs force per dir

close all;

% Compare Beam vs. Light Touch
% Plot vs. mean for t-test against mean
subplot(1,3,1)
boxplot(meanFx(:,2)); sigstar({1}); set(gca,'xticklabel',conds{4})
ylim([0 13]);
hline(1,'r-','1N');
box off; set(gca,'tickdir','out');
ylabel('Mean Fmag ML (N)');

subplot(1,3,2)
boxplot(meanFy(:,2)); sigstar({1}); set(gca,'xticklabel',conds{4})
ylim([0 13]);
hline(1,'r-','1N');
box off; set(gca,'tickdir','out');
ylabel('Mean Fmag AP (N)');

subplot(1,3,3)
boxplot(meanFz(:,2)); sigstar({1}); set(gca,'xticklabel',conds{4})
ylim([0 13]);
hline(1,'r-','1N');
box off; set(gca,'tickdir','out');
ylabel('Mean Fmag Vert (N)');

%% Plot mean abs force compare Overground vs. Beam
close all;
plotind = 0;

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanFx); %sigstar({[1,2]}); 
set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Fmag ML (N)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanFy); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Fmag AP (N)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanFz); sigstar({[1,2]}); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Fmag Vert (N)'); xtickangle(45)

%% Plot SD abs force compare Overground vs. Beam
close all;
plotind = 0;

boxplot(SDFx); sigstar({[1,2]}); 
set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('SD Fmag ML (N)'); xtickangle(45)

%%
close all
% Different color per subject
plotind = 0;

plotind = plotind + 1;
subplot(1,3,plotind)
plot(meanFx'); %sigstar({[1,2]}); 
set(gca,'xticklabel',conds([2 4]),'xtick',1:2); box off; set(gca,'tickdir','out');
ylabel('Mean Fmag ML (N)'); xtickangle(45), xlim([0.5 2.5]);
i = 0;
for subj = subj_array_force
    i = i + 1;
    legNames{i} = sprintf('HHI%i',subj);
end
legend(legNames);

plotind = plotind + 1;
subplot(1,3,plotind)
plot(meanFy'); set(gca,'xticklabel',conds([2 4]),'xtick',1:2); box off; set(gca,'tickdir','out');
ylabel('Mean Fmag AP (N)'); xtickangle(45), xlim([0.5 2.5]);

plotind = plotind + 1;
subplot(1,3,plotind)
plot(meanFz'); sigstar({[1,2]}); set(gca,'xticklabel',conds([2 4]),'xtick',1:2); box off; set(gca,'tickdir','out');
ylabel('Mean Fmag Vert (N)'); xtickangle(45), xlim([0.5 2.5]);

%% Plot mean abs IP vel per trial to compare with force mag's and power mag's
close all;
plotind = 0;

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanVx); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Vel Mag ML (m/s)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanVy); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
sigstar({[1 2]});
ylabel('Vel Mag AP (m/s)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanVz); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Vel Mag Vert. (m/s)'); xtickangle(45)

%% Plot mean abs, pos, and neg IP power overground vs. beam
plotind = 0;
numcols = 3; numrows = 1;

% % Abs power vs. thresh
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerAbsIntPtX(:,2)); sigstar({1}); 
% hline(0.0256,'k--'), ylim([0 0.21])
% set(gca,'xticklabel','Assist Beam'); box off; set(gca,'tickdir','out');
% ylabel('Power Mag ML (W)'); 
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerAbsIntPtY(:,2)); sigstar({1}); 
% hline(2.621,'k--'),ylim([0 2.5])
% set(gca,'xticklabel','Assist Beam'); box off; set(gca,'tickdir','out');
% ylabel('Power Mag AP (W)'); 
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerAbsIntPtX(:,2)); sigstar({1}); 
% hline(0.2597,'k--'),ylim([0 0.21]);
% set(gca,'xticklabel','Assist Beam'); box off; set(gca,'tickdir','out');
% ylabel('Power Mag Vert (W)'); 

%% Mean abs power (signif power flow)

plotind = 0;
plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(powerAbsIntPtX); %sigstar({[1,2]}); 
set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Power ML (W)'); xtickangle(45)

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(powerAbsIntPtY); sigstar({[1,2]})
set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Power AP (W)'); xtickangle(45)
title('Mean power magnitude Int. Pt.');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(powerAbsIntPtZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Power Vert (W)'); xtickangle(45)

%% Power mag variability

boxplot(SDpowerAbsIntPtX);
set(gca,'xticklabel',{'Partnered Ground','Partnered Beam'}); box off; set(gca,'tickdir','out');
ylabel('ML Power Variability (W)'); xtickangle(45)

%% Compare positive vs. negative power for Assist Beam in each dir (motor vs. brake)
figure
plotind = 0;

xlab{1} = 'Motor';
xlab{2} = 'Brake';
plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot([powerPosIntPtX(:,2) abs(powerNegIntPtX(:,2))]); 
% mPower = mean([powerPosIntPtX(:,2) powerNegIntPtX(:,2)]);
% sdPower = std([powerPosIntPtX(:,2) powerNegIntPtX(:,2)]);
% errorbar(1:2,mPower,sdPower)
sigstar({[1,2]});
set(gca,'xticklabel',xlab); box off; set(gca,'tickdir','out');
ylabel('ML Power (W)'); xtickangle(45)

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot([powerPosIntPtY(:,2) abs(powerNegIntPtY(:,2))]); sigstar({[1,2]});
set(gca,'xticklabel',xlab); box off; set(gca,'tickdir','out');
ylabel('Power AP (W)'); xtickangle(45)
title('Mean power Int. Pt. Assist Beam');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot([powerPosIntPtZ(:,2) abs(powerNegIntPtZ(:,2))]); sigstar({[1,2]});
set(gca,'xticklabel',xlab); box off; set(gca,'tickdir','out');
ylabel('Power Vert (W)'); xtickangle(45)

% % Positive power
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerPosIntPtX); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Mean power added Int Pt ML (W)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerPosIntPtY); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Mean power added Int Pt AP (W)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerPosIntPtZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Mean power added Int Pt Vert (W)'); xtickangle(45)
% 
% % Negative power
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(abs(powerNegIntPtX)); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Mean power absorb Int Pt ML (W)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(abs(powerNegIntPtY)); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Mean power absorb Int Pt AP (W)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(abs(powerNegIntPtZ)); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Mean power absorb Int Pt Vert (W)'); xtickangle(45)

%% Plot mean abs POB vel per trial to compare with force mag's and power mag's
% close all;
% plotind = 0;
% 
% plotind = plotind + 1;
% subplot(1,3,plotind)
% boxplot(meanVxPOB); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% sigstar({[1 2]});
% ylabel('Vel Mag ML (m/s)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(1,3,plotind)
% boxplot(meanVyPOB); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% sigstar({[1 2]});
% ylabel('Vel Mag AP (m/s)'); xtickangle(45)
% title('POB Clav');
% 
% plotind = plotind + 1;
% subplot(1,3,plotind)
% boxplot(meanVzPOB); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% sigstar({[1 2]});
% ylabel('Vel Mag Vert. (m/s)'); xtickangle(45)

%% Plot abs, pos, and neg power done on POB for trial for a bunch of conditions
% plotind = 0;
% numcols = 3; numrows = 1;
% 
% % Abs power (signif power flow)
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerAbsPOBX); 
% sigstar({[1,2]}); 
% set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Power ML (W)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerAbsPOBY); 
% sigstar({[1,2]})
% set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Power AP (W)'); xtickangle(45)
% title('Mean power magnitude POB');
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(powerAbsPOBZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Power Vert (W)'); xtickangle(45)
% 
% % Compare positive vs. negative power for Assist Beam in each dir
% figure
% plotind = 0;
% 
% xlab{1} = 'Added';
% xlab{2} = 'Absorbed';
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot([powerPosPOBX(:,2) abs(powerNegPOBX(:,2))]); 
% sigstar({[1,2]});
% set(gca,'xticklabel',xlab); box off; set(gca,'tickdir','out');
% ylabel('Power ML (W)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot([powerPosPOBY(:,2) abs(powerNegPOBY(:,2))]); 
% sigstar({[1,2]});
% set(gca,'xticklabel',xlab); box off; set(gca,'tickdir','out');
% ylabel('Power AP (W)'); xtickangle(45)
% title('Mean power POB Partnered Beam');
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot([powerPosPOBZ(:,2) abs(powerNegPOBZ(:,2))]); 
% sigstar({[1,2]});
% set(gca,'xticklabel',xlab); box off; set(gca,'tickdir','out');
% ylabel('Power Vert (W)'); xtickangle(45)

%% Concatenate performance metric data

clear; close all; clc;
load HHI2017_LateStats_MW.mat;
% Find means for each trial type for each subj
n = 0;
subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)
% conds_perf = {'Solo Beam','Assist Beam'};
conds_perf = {'Solo Ground','Assist Ground','Solo Beam','Assist Beam'};

n = 0;
for subj = subj_array % HHI12 had close to mean sway, good example subj
    n = n + 1;
    for i = 1:length(conds_perf)
        rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds_perf{i});
        StdSway(n,i) = nanmean(LateGroup.StdSway(rows)); % SD Clav sway
        Dist(n,i) = nanmean(LateGroup.Dist(rows)); % Total distance traveled on beam
        AvgSpeed(n,i) = nanmean(LateGroup.AvgSpeed(rows)); % Average speed of walking on beam
    end
end

%% Plot performance metric results, denote where significant diff.
close all;
subplot(1,3,1)
boxplot(StdSway); sigstar({[1,2]}); set(gca,'xticklabel',conds_perf); 
box off; set(gca,'tickdir','out');
ylabel('Std sway (mm)'); xtickangle(45)

% Plot vs. mean for t-test against mean
subplot(1,3,2)
boxplot(Dist(:,1)); sigstar({1}); set(gca,'xticklabel',conds_perf{1})
hline(Dist(1,2),'r-','beam length');
box off; set(gca,'tickdir','out');
ylabel('Distance traveled on beam (m)');

subplot(1,3,3)
boxplot(AvgSpeed); sigstar({[1,2]}); set(gca,'xticklabel',conds_perf); box off; set(gca,'tickdir','out');
ylabel('Avg. Speed (m/s)'); xtickangle(45)

%% Plot performance metric results, color code participants
close all;
subplot(1,3,1)
h = plot(StdSway'); sigstar({[1,2]});  
box off; set(gca,'tickdir','out'); xlim([0.5 2.5]); set(gca,'xtick',1:2);
set(gca,'xticklabel',conds_perf);
ylabel('Std sway (mm)'); xtickangle(45)

c = get(gcf,'colormap');

% Plot vs. mean for t-test against mean
subplot(1,3,2), hold on;
for n = 1:length(subj_array)
    plot(1,Dist(n,1),'x'); 
end
sigstar({1}); xlim([0.5 1.5]); set(gca,'xtick',1,'xticklabel',conds_perf{1})
hline(Dist(1,2),'r-','beam length');
box off; set(gca,'tickdir','out');
ylabel('Distance traveled on beam (m)');

subplot(1,3,3)
plot(AvgSpeed'); sigstar({[1,2]}); box off; set(gca,'tickdir','out');
xlim([0.5 2.5]); set(gca,'xtick',1:2);
set(gca,'xticklabel',conds_perf);
ylabel('Avg. Speed (m/s)'); xtickangle(45)

%% Plot kinematics for solo vs partner overground
close all;

n = 0;
subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)
conds_perf = {'Solo Ground','Assist Ground'};

n = 0;
for subj = subj_array % HHI12 had close to mean sway, good example subj
    n = n + 1;
    for i = 1:length(conds_perf)
        rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds_perf{i});
        StdSway(n,i) = nanmean(LateGroup.StdSway(rows)); % SD Clav sway
        AvgSpeed(n,i) = nanmean(LateGroup.AvgSpeed(rows)); % Average speed of walking on beam
    end
end

subplot(1,3,1)
boxplot(StdSway); set(gca,'xticklabel',conds_perf); 
box off; set(gca,'tickdir','out');
ylabel('ML sway (mm)'); xtickangle(45)

subplot(1,3,2)
boxplot(AvgSpeed); set(gca,'xticklabel',conds_perf); box off; set(gca,'tickdir','out');
ylabel('Mean forward speed (m/s)'); xtickangle(45)

%% Kinematics analysis: do for all subj's and cond's
for subj = subj_array(2:end-1) % only have 4-13 for angle data for now
    n = n + 1;
    % Do for all conditions (kinematic calculations)
    for i = 1:length(conds)
        rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds{i});
        % Frontal plane segment angles
%         StdPelvicObliq(n,i) = nanmean(LateGroup.StdPelvicObliq(rows));
        StdThoraxObliq(n,i) = nanmean(LateGroup.StdThoraxObliq(rows));
        StdLegObliq(n,i) = nanmean(LateGroup.StdLegObliq(rows));
%         rPelvicThoraxObliq(n,i) = nanmean(LateGroup.rPelvicThoraxObliq(rows));
        rLegThoraxObliq(n,i) = nanmean(LateGroup.rLegThoraxObliq(rows));
        % Sway metrics
        StdSway(n,i) = nanmean(LateGroup.StdSway(rows));
        StdCOMSway(n,i) = nanmean(LateGroup.StdCOMSway(rows));
        rFvertClav(n,i) = nanmean(LateGroup.rFvertClav(rows));
        rFlatCOM(n,i) = nanmean(LateGroup.rFlatCOM(rows));
        rFvertCOM(n,i) = nanmean(LateGroup.rFvertCOM(rows));
        rFlatClav(n,i) = nanmean(LateGroup.rFlatClav(rows));
    end
end

%% Plot leg and thorax obliquity angles

subplot(2,2,1)
% boxplot(StdPelvicObliq); set(gca,'xticklabel',conds); box off; set(gca,'tickdir','out');
% ylabel('Std pelvic obliq (deg)'); xtickangle(45)
boxplot(StdLegObliq); set(gca,'xticklabel',conds); box off; set(gca,'tickdir','out');
ylabel('Std Leg obliq (deg)'); xtickangle(45)

subplot(2,2,2)
boxplot(StdThoraxObliq); set(gca,'xticklabel',conds); box off; set(gca,'tickdir','out');
ylabel('Std thorax obliq (deg)'); xtickangle(45)

subplot(2,2,3)
boxplot(rLegThoraxObliq); set(gca,'xticklabel',conds); box off; set(gca,'tickdir','out');
ylabel('Corr Leg/thorax obliq'); xtickangle(45)

% subplot(2,2,3)
% boxplot(rPelvicThoraxObliq); set(gca,'xticklabel',conds); box off; set(gca,'tickdir','out');
% ylabel('Corr pelvic/thorax obliq'); xtickangle(45)

%% First compare if Std Clav sway is very different from Std COM sway

subplot(1,2,1)
boxplot(StdSway); set(gca,'xticklabel',conds); box off; set(gca,'tickdir','out');
ylabel('Std sway (mm)'); xtickangle(45)

subplot(1,2,2)
boxplot(StdCOMSway); set(gca,'xticklabel',conds); box off; set(gca,'tickdir','out');
ylabel('Std COM sway (mm)'); xtickangle(45)

%% Compare if correspondence between force and sway differ depending on which type of F or sway data used
figure;
assistConds = {'Assist Ground';'Assist Beam'};
subplot(2,2,1)
boxplot(rFvertClav(:,[2 4])); set(gca,'xticklabel',assistConds); box off; set(gca,'tickdir','out');
ylabel('Corr. Clav sway to vert. F'), ylim([-1 1]);

subplot(2,2,2)
boxplot(rFlatCOM(:,[2 4])); set(gca,'xticklabel',assistConds); box off; set(gca,'tickdir','out');
ylabel('Corr. COM sway to lat. F'), ylim([-1 1]);

subplot(2,2,3)
boxplot(rFvertCOM(:,[2 4])); set(gca,'xticklabel',assistConds); box off; set(gca,'tickdir','out');
ylabel('Corr. COM sway to vert. F'), ylim([-1 1]);

subplot(2,2,4)
boxplot(rFlatClav(:,[2 4])); set(gca,'xticklabel',assistConds); box off; set(gca,'tickdir','out');
ylabel('Corr. Clav sway to lat. F'), ylim([-1 1]);

%% Old plots after this point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot mean net work done per trial

plotind = 0;

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanWorkIntPtX); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Work ML (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanWorkIntPtY); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Work AP (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanWorkIntPtZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Work Vert (J)'); xtickangle(45)

%% Plot mean abs work done per trial - boxplots

plotind = 0;

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanWorkMagIntPtX); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Abs Work ML (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanWorkMagIntPtY); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Abs Work AP (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanWorkMagIntPtZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Abs Work Vert (J)'); xtickangle(45)

%% Plot mean abs work done per trial - color code participants

plotind = 0;

plotind = plotind + 1;
subplot(1,3,plotind)
plot(meanWorkMagIntPtX'); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Abs Work ML (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanWorkMagIntPtY); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Abs Work AP (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(1,3,plotind)
boxplot(meanWorkMagIntPtZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Abs Work Vert (J)'); xtickangle(45)

%% Plot mean pos and neg work done per trial

plotind = 0;

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosWorkIntPtX); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos Work ML (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosWorkIntPtY); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos Work AP (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosWorkIntPtZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos Work Vert (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegWorkIntPtX); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg Work ML (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegWorkIntPtY); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg Work AP (J)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegWorkIntPtZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg Work Vert (J)'); xtickangle(45)

%% Plot mean pos and neg work done by assistant at interaction pt
% All forces are resolved into lab coordinate frame, tension is positive and in neg z dir w.r.t.
% Assistant's finger/arm. If force x disp. is negative in lab CS, positive
% work done by assistant.

% subplot(2,2,1)
% boxplot(PosWorkAssist); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Peak pos. work Asst. arm len (J)'); xtickangle(45)
% 
% subplot(2,2,2)
% boxplot(NegWorkAssist); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Peak neg. work Asst. arm len (J)'); xtickangle(45)

%% Plot total positive and negative work done for condition - biased by time length of trial
% plotind = 0;
% 
% plotind = plotind + 1;
% subplot(2,4,plotind)
% boxplot(PosWorkIntPtTot); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Energy added Int Pt Vec (J)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(2,4,plotind)
% boxplot(PosWorkIntPtX); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Energy added Int Pt ML (J)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(2,4,plotind)
% boxplot(PosWorkIntPtY); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Energy added Int Pt AP (J)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(2,4,plotind)
% boxplot(PosWorkIntPtZ); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Energy added Int Pt Vert (J)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(2,4,plotind)
% boxplot(abs(NegWorkIntPtTot)); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Energy absorbed Int Pt Vec (J)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(2,4,plotind)
% boxplot(abs(NegWorkIntPtX)); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Energy absorbed Int Pt ML (J)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(2,4,plotind)
% boxplot(abs(NegWorkIntPtY)); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Energy absorbed Int Pt AP (J)'); xtickangle(45)
% 
% plotind = plotind + 1;
% subplot(2,4,plotind)
% boxplot(abs(NegWorkIntPtZ)); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
% ylabel('Energy absorbed Int Pt Vert (J)'); xtickangle(45)

 %% Plot regression IP R^2 and coeff values ML dir Partner Beam only
% plotind = 0;
% 
% plotind = plotind + 1;
% subplot(4,1,plotind)
% boxplot(meanRsqx(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Rsq'); title('ML direction');
% 
% plotind = plotind + 1;
% subplot(4,1,plotind)
% boxplot(meanMx(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Mass (kg)'); 
% 
% plotind = plotind + 1;
% subplot(4,1,plotind)
% boxplot(meanBx(:,2)); set(gca,'xticklabel','All partnerships');box off; set(gca,'tickdir','out');
% ylabel('Damping (Ns/m)'); 
% 
% plotind = plotind + 1;
% subplot(4,1,plotind)
% boxplot(meanKx(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Stiffness (N/m)'); 
% 
% %% Plot regression IP R^2 and coeff values Vert dir Partner Beam only
% plotind = 0;
% 
% plotind = plotind + 1;
% subplot(4,1,plotind)
% boxplot(meanRsqz(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('R^2'); title('Vertical direction');
% 
% plotind = plotind + 1;
% subplot(4,1,plotind)
% boxplot(meanMz(:,2)); hline(0,'k--'),
% set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('m (kg)'); 
% 
% plotind = plotind + 1;
% subplot(4,1,plotind)
% boxplot(meanBz(:,2)), hline(0,'k--'); 
% set(gca,'xticklabel','All partnerships');box off; set(gca,'tickdir','out');
% ylabel('b (Ns/m)'); 
% 
% plotind = plotind + 1;
% subplot(4,1,plotind)
% boxplot(meanKz(:,2)),hline(0,'k--');
% set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('k (N/m)'); 

%% mean percent Partner Beam trials where m, b, k sig. ML dir IP

plotind = 0;
xlab = 'All partnerships';
numrows = 3; numcols = 1;

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(mxIPsig*100),ylabel('% trials m sig')
box off; ylim([0 100]);
set(gca,'xticklabel',xlab,'tickdir','out');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(bxIPsig*100),ylabel('% trials b sig')
box off; ylim([0 100]);
set(gca,'xticklabel',xlab,'tickdir','out');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(kxIPsig*100),ylabel('% trials k sig')
box off; ylim([0 100]);
set(gca,'xticklabel',xlab,'tickdir','out');

%% Plot mean pos and neg force
close all;
plotind = 0;

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosFx); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos F ML (N)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosFy); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos F AP (N)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosFz); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos F Vert (N)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegFx); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg F ML (N)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegFy); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg F AP (N)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegFz); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg F Vert (N)'); xtickangle(45)

%% Plot mean pos and neg IP vel per trial
figure;
plotind = 0;

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosVx); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos Vel ML (m/s)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosVy); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos Vel AP (m/s)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanPosVz); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Pos Vel Vert (m/s)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegVx); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg Vel ML (m/s)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegVy); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg Vel AP (m/s)'); xtickangle(45)

plotind = plotind + 1;
subplot(2,3,plotind)
boxplot(meanNegVz); set(gca,'xticklabel',conds([2 4])); box off; set(gca,'tickdir','out');
ylabel('Mean Neg Vel Vert (m/s)'); xtickangle(45)
