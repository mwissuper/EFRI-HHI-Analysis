%% Plot results for all subj's combined
clear; close all; clc;

colors = hsv(12);
colors(3,:) = [0.85,0.33,0.10]; % replace yellow
colors(5,:) = [0.47,0.67,0.19]; % replace a green
colors(6,:) = [0.5,0.5,0.5]; 
colors(8,:) = [0.00,0.50,1.00];
colors(12,:) = [0.49,0.18,0.56];

%% Calculate subj means for performance metric data

load HHI2017_Stats_MW.mat;
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
        rows = GroupStats.Subject==subj & strcmp(GroupStats.Type,conds_perf{i});
        StdSway(n,i) = nanmean(GroupStats.StdSway(rows)); % SD Clav sway
        Dist(n,i) = nanmean(GroupStats.Dist(rows)); % Total distance traveled on beam
        AvgSpeed(n,i) = nanmean(GroupStats.AvgSpeed(rows)); % Average speed of walking on beam
        
        meanLyGd(n,i) = nanmean(GroupStats.meanLy(rows)); % Just mean, not mean of abs value
        SDLyGd(n,i) = nanmean(GroupStats.SDLy(rows));
        
        meanW(n,i) = nanmean(GroupStats.meanW(rows)); % Just mean, not mean of abs value
    end
end

%% Plot performance metric results, color code participants
close all;

xlab = {'Solo ground','Partnered ground','Solo beam','Partnered beam'};

% Distance on beam
% Plot vs. mean for t-test against mean
subplot(1,3,1)
hold on;

% % Plot each subj's data
% for n = 1:length(subj_array)
%     if n == 2 
%         plot(-.25,Dist(n,1),'.','markersize',14,'color',colors(n,:)); 
%     elseif n == 4 % similar value to n = 2
%         plot(.25,Dist(n,1),'.','markersize',14,'color',colors(n,:)); 
%     elseif n == 10 
%         plot(-.25,Dist(n,1),'.','markersize',14,'color',colors(n,:)); 
%     elseif n == 12 
%         plot(.25,Dist(n,1),'.','markersize',14,'color',colors(n,:));  
%     elseif n == 6 
%         plot(-.25,Dist(n,1),'.','markersize',14,'color',colors(n,:)); 
%     elseif n == 7 
%         plot(-.25,Dist(n,1),'.','markersize',14,'color',colors(n,:)); 
%     elseif n == 8 
%         plot(.25,Dist(n,1),'.','markersize',14,'color',colors(n,:)); 
%     elseif n == 11 % similar value to n = 10 and 12, space apart on x axis
%         plot(0,Dist(n,1),'.','markersize',14,'color',colors(n,:)); 
%     else
%         plot(0,Dist(n,1),'.','markersize',14,'color',colors(n,:)); 
%     end
% end
% boxplot(Dist(:,1),zeros(size(Dist(:,1))),'colors','k','symbol','o','widths',0.25); %sigstar({1}); set(gca,'xticklabel',conds_perf{1})
% xlim([0 2]); set(gca,'xtick',5)
bar(1:2,[mean(Dist(:,3)) 3.65],'linestyle','none')
errorbar(1,mean(Dist(:,3)),nanstd(Dist(:,3)),'Linestyle','none','color','k');
% axis tight
ylim([0 3.7]),xlim([0.5 2.5]),set(gca,'xtick',1:4,'xticklabel',{'Solo beam','Partnered beam'})
hline(3.65,'k--','beam length');
box off; set(gca,'tickdir','out');
% title('Solo distance completed'),
title('Distance completed'),ylabel('(m)');

% Sway
subplot(1,3,2),hold on
% for n = 1:length(subj_array)
%     plot(1:2,StdSway(n,:),'color',colors(n,:)); 
% end
bar([1 2 4 5],mean(StdSway),'linestyle','none')
errorbar([1 2 4 5],mean(StdSway),nanstd(StdSway),'Linestyle','none','color','k');
sigstar({[4,5]});  
% boxplot(StdSway,'colors','k','symbol','o'); %sigstar({[1,2]}); set(gca,'xticklabel',condLab); 
box off; set(gca,'tickdir','out'); 
% xlim([0.5 2.5]); set(gca,'xtick',1:2);
% set(gca,'xticklabel',{'Solo','Partnered'}); 
xlim([0.5 5.5]),ylim([0 0.12]),set(gca,'xtick',[1 2 4 5],'xticklabel',xlab)
set(gca,'ytick',0:0.04:0.12)
title('Sway variability'),ylabel('(m)'); 

% Plot change in metric vs solo beam distance
% subplot(2,3,4),hold on
% c = polyfit(Dist(:,1),AvgSpeed(:,2) - AvgSpeed(:,1),1);
% for n = 1:length(subj_array)
%     plot(Dist(n,1),AvgSpeed(n,2) - AvgSpeed(n,1),'.','markersize',14,'color',colors(n,:)); 
% end
% % plot(Dist(:,1),polyval(c,Dist(:,1)),'k');
% xlabel('Solo distance completed');  %ylim([-0.2 .5])
% box off; set(gca,'tickdir','out');
% title('Change in forward speed'),ylabel('(m/s)');
[p,r,test] = corr2vars(Dist(:,3),StdSway(:,3)-StdSway(:,4))
% Change in sway var vs. solo beam performance
subplot(1,3,3),hold on;
for n = 1:length(subj_array)
    plot(Dist(n,3),StdSway(n,3)-StdSway(n,4),'.','markersize',14,'color',colors(n,:)); 
end
xlabel('Solo distance completed'); 
box off; set(gca,'tickdir','out');xlim([0 3.7])
title('Reduction in sway variability'),ylabel('(m)');ylim([0 0.085])
set(gca,'ytick',0:0.04:0.08)

% % Plot change in metric vs solo beam distance
% subplot(2,3,4)
% c = polyfit(Dist(:,1),AvgSpeed(:,2) - AvgSpeed(:,1),1);
% plot(Dist(:,1),AvgSpeed(:,2) - AvgSpeed(:,1),'kx'),hold on;
% plot(Dist(:,1),polyval(c,Dist(:,1)),'k');
% xlabel('Solo distance completed'); 
% box off; set(gca,'tickdir','out');
% ylabel('Change in forward speed (m/s)');
% 
% subplot(2,3,5)
% plot(Dist(:,1),StdSway(:,2)-StdSway(:,1),'kx'),
% xlabel('Solo distance completed'); 
% box off; set(gca,'tickdir','out');
% ylabel('Change in sway (m/s)'); 

% set(gcf,'outerposition',[267   423   685   457]);

%% Plot mean angular momentum and correlation with sway metric
% bar plot
figure
numrows = 1; numcols = 2;
plotind = 0;
condLab{1} = 'Solo';
condLab{2} = 'Partnered';

% Means and SD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean Ly about ground axis
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanLyGd(:,3:4),condLab,'Mean RMS Ly gd','(kg*m^2/s)');
sigstar({[1,2]})

% SD Ly about ground axis
% subplot(1,2,2),hold on;
% h = barPlot2cond(SDLyGd(:,3:4),condLab,'Ly gd variability','(kg*m/s)');

clear p r test ind

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
ind = 1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr([meanLyGd(:,3); meanLyGd(:,4)],[StdSway(:,3);StdSway(:,4)],subj_array,subj_array,colors,'Sway vs. Ly gd','(rad/s)','SD ML torso disp. (m)'); 

%% Do same as above for angular speed (RMS) just to double check

% bar plot
figure
numrows = 1; numcols = 2;
plotind = 0;
condLab{1} = 'Solo';
condLab{2} = 'Partnered';

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanW(:,3:4),condLab,'Mean RMS torso ang speed','(rad/s)');

clear p r test ind

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
ind = 1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr([meanW(:,3); meanW(:,4)],[StdSway(:,3);StdSway(:,4)],subj_array,subj_array,colors,'Sway vs. W','(rad/s)','SD ML torso disp. (m)'); 

%% Plot AP speed (strategy)
% Speed
hold on
for n = 1:length(subj_array)
    plot(1:2,AvgSpeed(n,3:4),'color',colors(n,:)); 
end
box off; set(gca,'tickdir','out');
xlim([0.5 2.5]); set(gca,'xtick',1:2);
set(gca,'xticklabel',{'Solo','Partnered'});
title('Mean forward speed');ylabel('(m/s)')

boxplot(AvgSpeed,'colors','k','symbol','o'); % sigstar({[1,2]}); set(gca,'xticklabel',condLab); box off; set(gca,'tickdir','out');
% ylabel('Mean forward speed (m/s)'); xtickangle(45)
xlim([0.5 2.5]); set(gca,'xtick',1:2);box off;
set(gca,'xticklabel',{'Solo','Partnered'}); %ylim([-0.05 0.8]),
% set(gca,'ytick',0:0.4:0.8)

%% Force data

load HHI2017_Stats_force_MW.mat;
% Find means for each trial type for each subj
n = 0;

subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)
colors = hsv(12);
colors(3,:) = [0.85,0.33,0.10]; % replace yellow
colors(5,:) = [0.47,0.67,0.19]; % replace a green
colors(6,:) = [0.5,0.5,0.5]; 
colors(8,:) = [0.00,0.50,1.00];
colors(12,:) = [0.49,0.18,0.56];
subj_array = 3:14;
conds = {'Solo Ground','Assist Ground','Solo Beam','Assist Beam'};

%% Get subset of performance data that corrseponds to subj's with force data

[c,ia,ib] = intersect(subj_array,subj_array_force);
soloDistF = Dist(ia,1); % solo ability
SwayVRedF = StdSway(ia,3) - StdSway(ia,4); % improvement in balance

%% Calculate subj means for force, vel, power data: do only for subj's with good force data and for assist cond's

n = 0;
for subj = subj_array_force
    n = n + 1;
    for i = [2 4]
        rows = GroupStats.Subject==subj & strcmp(GroupStats.Type,conds{i});
        
        %% ML/X dir
        
        % Force 
        meanFx(n,i/2) = nanmean(GroupStats.meanFx(rows)); % Just mean, not mean of abs value
        SDFx(n,i/2) = nanmean(GroupStats.SDFx(rows));
        % Velocity of IP
        meanVx(n,i/2) = nanmean(GroupStats.meanVx(rows)); % Mean abs(v)
        SDVx(n,i/2) = nanmean(GroupStats.SDVx(rows)); % Mean abs(v)
        % Power at IP 
        meanPX(n,i/2) = nanmean(GroupStats.meanIPpowerX(rows)); % signed mean
        SDPX(n,i/2) = nanmean(GroupStats.SDIPpowerX(rows));
        % Regression beam-walker torso to force (nan value if n.s. fit)
        meanMxPOB(n,i/2) = nanmean(GroupStats.mx_torso(rows));
        meanBxPOB(n,i/2) = nanmean(GroupStats.bx_torso(rows));
        meanKxPOB(n,i/2) = nanmean(GroupStats.kx_torso(rows));
        meanRsqxPOB(n,i/2) = nanmean(GroupStats.Rsqx_torso(rows));
        % Calculate percent of trials where param was sig.
        if i == 4 % Assist Beam condition only
            mxPOBsig(n) = length(find(~isnan(GroupStats.mx_torso(rows))))/length(GroupStats.mx_torso(rows));
            bxPOBsig(n) = length(find(~isnan(GroupStats.bx_torso(rows))))/length(GroupStats.bx_torso(rows));
            kxPOBsig(n) = length(find(~isnan(GroupStats.kx_torso(rows))))/length(GroupStats.kx_torso(rows));
        end
        
        %% Vertical/Z dir
        
        % Force 
        meanFz(n,i/2) = nanmean(GroupStats.meanFz(rows)); % Just mean, not mean of abs value
        SDFz(n,i/2) = nanmean(GroupStats.SDFz(rows));
        % Velocity of IP
        meanVz(n,i/2) = nanmean(GroupStats.meanVz(rows)); % Mean abs(v)
        SDVz(n,i/2) = nanmean(GroupStats.SDVz(rows)); % Mean abs(v)
        % Power at IP 
        meanPZ(n,i/2) = nanmean(GroupStats.meanIPpowerZ(rows)); % signed mean
        SDPZ(n,i/2) = nanmean(GroupStats.SDIPpowerZ(rows));
        % Regression beam-walker torso to force (nan value if n.s. fit)
        meanMzPOB(n,i/2) = nanmean(GroupStats.mz_torso(rows));
        meanBzPOB(n,i/2) = nanmean(GroupStats.bz_torso(rows));
        meanKzPOB(n,i/2) = nanmean(GroupStats.kz_torso(rows));
        meanRsqzPOB(n,i/2) = nanmean(GroupStats.Rsqz_torso(rows));
        % Calculate percent of trials where param was sig.
        if i == 4 % Assist Beam condition only
            mzPOBsig(n) = length(find(~isnan(GroupStats.mz_torso(rows))))/length(GroupStats.mz_torso(rows));
            bzPOBsig(n) = length(find(~isnan(GroupStats.bz_torso(rows))))/length(GroupStats.bz_torso(rows));
            kzPOBsig(n) = length(find(~isnan(GroupStats.kz_torso(rows))))/length(GroupStats.kz_torso(rows));
        end

        %% 2D vector force in x and z
        % Force 
        meanFxz(n,i/2) = nanmean(GroupStats.meanFxz(rows)); % Just mean, not mean of abs value
        SDFxz(n,i/2) = nanmean(GroupStats.SDFxz(rows));
        meanFxz(n,i/2) = nanmean(GroupStats.meanFxz(rows)); % Just mean, not mean of abs value
        SDFxz(n,i/2) = nanmean(GroupStats.SDFxz(rows));
        meanTheta(n,i/2) = nanmean(GroupStats.meanTheta(rows)).*180/pi; % Convert to deg
        SDtheta(n,i/2) = nanmean(GroupStats.SDTheta(rows)).*180/pi; % Convert to deg
%         % Velocity of IP
%         meanVxz(n,i/2) = nanmean(GroupStats.meanVxz(rows)); % Mean abs(v)
%         SDVxz(n,i/2) = nanmean(GroupStats.SDVxz(rows)); % Mean abs(v)
        % Power at IP 
        meanPXZ(n,i/2) = nanmean(GroupStats.meanIPpowerXZ(rows)); % signed mean
        SDPXZ(n,i/2) = nanmean(GroupStats.SDIPpowerXZ(rows));
        
        %% Torque
        meanTyGd(n,i/2) = nanmean(GroupStats.meanTorqueYgd(rows)); % RMS mean
        SDTyGd(n,i/2) = nanmean(GroupStats.SDTorqueYgd(rows));
        meanTyTorso(n,i/2) = nanmean(GroupStats.meanTorqueYTorso(rows)); % RMS mean
        SDTyTorso(n,i/2) = nanmean(GroupStats.SDTorqueYTorso(rows));
        
        %% Angular velocity and power
        meanTorsoW(n,i/2) = nanmean(GroupStats.meanW(rows)); % RMS mean
        meanPangGd(n,i/2) = nanmean(GroupStats.meanPowerAng(rows)); % signed mean
        meanPangGdRMS(n,i/2) = nanmean(GroupStats.meanPowerAngRMS(rows)); % RMS mean
        meanPangGdPos(n,i/2) = nanmean(GroupStats.meanPangPos(rows));
        meanPangGdNeg(n,i/2) = nanmean(GroupStats.meanPangNeg(rows));
        
        %% Torque cross-corr to torso state
        meanLagTyAngTorso(n,i/2) = nanmean(GroupStats.lagTyAngTorso(rows));
        meanXcorrTyAngTorso(n,i/2) = nanmean(GroupStats.xcorrTyAngTorso(rows));
        meanLagTyWTorso(n,i/2) = nanmean(GroupStats.lagTyWTorso(rows));
        meanXcorrTyWTorso(n,i/2) = nanmean(GroupStats.xcorrTyWTorso(rows));
        meanLagTyAlphaTorso(n,i/2) = nanmean(GroupStats.lagTyAlphaTorso(rows));
        meanXcorrTyAlphaTorso(n,i/2) = nanmean(GroupStats.xcorrTyAlphaTorso(rows));
        
        %% Model fit rotational
        % Regression beam-walker torso to force (nan value if n.s. fit)
        meanMangPOB(n,i/2) = nanmean(GroupStats.m_ang_torso(rows));
        meanBangPOB(n,i/2) = nanmean(GroupStats.b_ang_torso(rows));
        meanKangPOB(n,i/2) = nanmean(GroupStats.k_ang_torso(rows));
        meanRsqangPOB(n,i/2) = nanmean(GroupStats.Rsq_ang_torso(rows));
        % Calculate percent of trials where param was sig.
        if i == 4 % Assist Beam condition only
            mangPOBsig(n) = length(find(~isnan(GroupStats.m_ang_torso(rows))))/length(GroupStats.m_ang_torso(rows));
            bangPOBsig(n) = length(find(~isnan(GroupStats.b_ang_torso(rows))))/length(GroupStats.b_ang_torso(rows));
            kangPOBsig(n) = length(find(~isnan(GroupStats.k_ang_torso(rows))))/length(GroupStats.k_ang_torso(rows));
        end
        
        %% Old metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % corr power to Fresid
        meanrPowerFresX(n,i/2) = nanmean(GroupStats.rPowerFresX(rows));
        % Pos and neg vel
        meanPosVx(n,i/2) = nanmean(GroupStats.meanPosVx(rows));
        meanNegVx(n,i/2) = nanmean(GroupStats.meanNegVx(rows));
        % Pos and neg power and force
        meanPosFx(n,i/2) = nanmean(GroupStats.meanPosFx(rows));
        meanNegFx(n,i/2) = nanmean(GroupStats.meanNegFx(rows));
        SDPosFx(n,i/2) = nanmean(GroupStats.SDPosFx(rows));
        SDNegFx(n,i/2) = nanmean(GroupStats.SDNegFx(rows));
        meanPposX(n,i/2) = nanmean(GroupStats.meanPosPowerIntPtX(rows));
        meanPnegX(n,i/2) = nanmean(GroupStats.meanNegPowerIntPtX(rows));
        perPpos(n,i/2) = nanmean(GroupStats.perPposFposX(rows) + GroupStats.perPposFnegX(rows));
        perPneg(n,i/2) = nanmean(GroupStats.perPnegFposX(rows) + GroupStats.perPnegFnegX(rows));
        perVpos(n,i/2) = nanmean(GroupStats.perVposFposX(rows) + GroupStats.perVposFnegX(rows));
        perVneg(n,i/2) = nanmean(GroupStats.perVnegFposX(rows) + GroupStats.perVnegFnegX(rows));
        perVposFlo(n,i/2) = nanmean(GroupStats.perVposFloposX(rows) + GroupStats.perVposFlonegX(rows));
        perVnegFlo(n,i/2) = nanmean(GroupStats.perVnegFloposX(rows) + GroupStats.perVnegFlonegX(rows));
        perVposFhi(n,i/2) = nanmean(GroupStats.perVposFhiposX(rows) + GroupStats.perVposFhinegX(rows));
        perVnegFhi(n,i/2) = nanmean(GroupStats.perVnegFhiposX(rows) + GroupStats.perVnegFhinegX(rows));
        meanPposF(n,i/2) = nanmean(GroupStats.meanIPpowerPosF(rows));
        meanPnegF(n,i/2) = nanmean(GroupStats.meanIPpowerNegF(rows));
        meanPposFlo(n,i/2) = nanmean(GroupStats.meanIPpowerPosFlo(rows));
        meanPnegFlo(n,i/2) = nanmean(GroupStats.meanIPpowerNegFlo(rows));
        meanPposFhi(n,i/2) = nanmean(GroupStats.meanIPpowerPosFhi(rows));
        meanPnegFhi(n,i/2) = nanmean(GroupStats.meanIPpowerNegFhi(rows));
        % xcorr IP to POB Torso
        meanxcorrFIPTorsoX(n,i/2) = nanmean(GroupStats.xcorrFIPTorsoX(rows));
        meanxcorrFIPvTorsoX(n,i/2) = nanmean(GroupStats.xcorrFIPvTorsoX(rows));
        meanxcorrvIPvTorsoX(n,i/2) = nanmean(GroupStats.xcorrvIPvTorsoX(rows));
        % Force strategy proportion of time
        perFpos(n,i/2) = nanmean(GroupStats.perPposFposX(rows) + GroupStats.perPnegFposX(rows));
        perFneg(n,i/2) = nanmean(GroupStats.perPposFnegX(rows) + GroupStats.perPnegFnegX(rows));
        % xcorr lag IP to POB Torso
        meanLagFIPTorsoX(n,i/2) = nanmean(GroupStats.lagFIPTorsoX(rows));
        meanLagFIPvTorsoX(n,i/2) = nanmean(GroupStats.lagFIPvTorsoX(rows));
        meanLagvIPvTorsoX(n,i/2) = nanmean(GroupStats.lagvIPvTorsoX(rows));
        % POB armLen (not force but related to xcorr above)
        meanArmPOBX(n,i/2) = nanmean(GroupStats.meanArmPOBX(rows));
        SDarmPOBX(n,i/2) = nanmean(GroupStats.SDarmPOBX(rows));
        
        % Calculate corr's between pos/neg power and b/k values using all 
        % trials per subj. Only do for assist beam trials. Put each corr 
        % pair in column of metric. Col's: 1 = pos power and b, 2 = pos
        % power and k, 3 = neg power and b, 4 = neg power and k. corr done
        % within each subj. Later code uses all subj and trials together
        if i == 4
            [r,p] = corr(GroupStats.meanPosPowerIntPtX(rows),GroupStats.bx_torso(rows));
            if p < 0.05
                corrPbk(n,1) = r;
            else
                corrPbk(n,1) = nan;
            end
            
            [r,p] = corr(GroupStats.meanPosPowerIntPtX(rows),GroupStats.kx_torso(rows));
            if p < 0.05
                corrPbk(n,2) = r;
            else
                corrPbk(n,2) = nan;
            end
            
            [r,p] = corr(GroupStats.meanNegPowerIntPtX(rows),GroupStats.bx_torso(rows));
            if p < 0.05
                corrPbk(n,3) = r;
            else
                corrPbk(n,3) = nan;
            end
            
            [r,p] = corr(GroupStats.meanNegPowerIntPtX(rows),GroupStats.kx_torso(rows));
            if p < 0.05
                corrPbk(n,4) = r;
            else
                corrPbk(n,4) = nan;
            end
        end

        % Model coefficients POB torso x/ML dir
        % Keep all values for all trials to do boxplot for each subj
        if i == 4 % do for assist beam only
            meanMxArray{n} = GroupStats.mx_torso(rows);
            meanBxArray{n} = GroupStats.bx_torso(rows);
            meanKxArray{n} = GroupStats.kx_torso(rows);
            meanRsqxArray{n} = GroupStats.Rsqx_torso(rows);
            % with lag
            meanMxLagArray{n} = GroupStats.mx_lag_torso(rows);
            meanBxLagArray{n} = GroupStats.bx_lag_torso(rows);
            meanKxLagArray{n} = GroupStats.kx_lag_torso(rows);
            meanRsqxLagArray{n} = GroupStats.Rsqx_lag_torso(rows);
        end
        
        % torso x/ML dir
        % Velocity torso marker
        meanVxPOB(n,i/2) = nanmean(GroupStats.meanVx_torso(rows));
        % Standardized coefficients
        % Regression beam-walker torso to force (nan value if n.s. fit)
        meanMxPOBst(n,i/2) = nanmean(GroupStats.mx_torso_st(rows));
        meanBxPOBst(n,i/2) = nanmean(GroupStats.bx_torso_st(rows));
        meanKxPOBst(n,i/2) = nanmean(GroupStats.kx_torso_st(rows));
        meanRsqxPOBst(n,i/2) = nanmean(GroupStats.Rsqx_torso_st(rows));
        % lag
        meanMxLagPOB(n,i/2) = nanmean(GroupStats.mx_lag_torso(rows));
        meanBxLagPOB(n,i/2) = nanmean(GroupStats.bx_lag_torso(rows));
        meanKxLagPOB(n,i/2) = nanmean(GroupStats.kx_lag_torso(rows));
        meanRsqxLagPOB(n,i/2) = nanmean(GroupStats.Rsqx_lag_torso(rows));

        % Regression rFin to force (nan value if n.s. fit) ML dir
        meanMxIP(n,i/2) = nanmean(GroupStats.mx_IP(rows));
        meanBxIP(n,i/2) = nanmean(GroupStats.bx_IP(rows));
        meanKxIP(n,i/2) = nanmean(GroupStats.kx_IP(rows));
        meanRsqxIP(n,i/2) = nanmean(GroupStats.Rsqx_IP(rows));
        
        % Calculate subj mean and SD P and F sign combo, Assist Beam cond only
        if i == 4
            perPFmean(n,1) = nanmean(GroupStats.perPposFposX(rows)); % Q1
            perPFmean(n,4) = nanmean(GroupStats.perPposFnegX(rows)); % Q4
            perPFmean(n,2) = nanmean(GroupStats.perPnegFposX(rows)); % Q2
            perPFmean(n,3) = nanmean(GroupStats.perPnegFnegX(rows)); % Q3
            perPFSD(n,1) = nanstd(GroupStats.perPposFposX(rows)); % Q1
            perPFSD(n,4) = nanstd(GroupStats.perPposFnegX(rows)); % Q4
            perPFSD(n,2) = nanstd(GroupStats.perPnegFposX(rows)); % Q2
            perPFSD(n,3) = nanstd(GroupStats.perPnegFnegX(rows)); % Q3
        end
    end
end

save LateStats_groupMeans_force

%% Fx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot means and SD's for IP force, vel, and power Overground vs. Beam x dir
% bar plots
% figure
condLab{1} = 'Ground';
condLab{2} = 'Beam';

% Means and SD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean Force IP
subplot(2,6,1),hold on;
bar(mean(meanFx),'linestyle','none');
errorbar(mean(meanFx),std(meanFx),'linestyle','none','color','k');
% boxplot(meanFx,'colors','k','symbol','o')
set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('Mean force'),ylabel('(N)'); %xtickangle(45)
xlim([0.5 2.5]);

% SD Force IP
subplot(2,6,2),hold on;
bar(mean(SDFx),'linestyle','none');
errorbar(mean(SDFx),std(SDFx),'linestyle','none','color','k');
% boxplot(SDFx,'colors','k','symbol','o')
sigstar({[1,2]});
set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('Force variability'),ylabel('(N)'); %xtickangle(45
xlim([0.5 2.5]);

% Mean Vel IP
subplot(2,6,3),hold on;
bar(mean(meanVx),'linestyle','none');
errorbar(mean(meanVx),std(meanVx),'linestyle','none','color','k');
% boxplot(meanVx,'colors','k','symbol','o')
set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('Mean RMS velocity'),ylabel('(m/s)'); %xtickangle(45)
xlim([0.5 2.5]);

% SD Vel IP
subplot(2,6,4),hold on;
bar(mean(SDVx),'linestyle','none');
errorbar(mean(SDVx),std(SDVx),'linestyle','none','color','k');
% boxplot(SDVx,'colors','k','symbol','o')
set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('RMS velocity variability'),ylabel('(m/s)'); %xtickangle(45
xlim([0.5 2.5]);

% Mean Power IP
subplot(2,6,5),hold on;
bar(mean(meanPX),'linestyle','none');
errorbar(mean(meanPX),std(meanPX),'linestyle','none','color','k');
% boxplot(meanPX,'colors','k','symbol','o')
set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('Mean power'),ylabel('(W)'); %xtickangle(45
xlim([0.5 2.5]);

% SD Power IP
subplot(2,6,6),hold on;
bar(mean(SDPX),'linestyle','none');
errorbar(mean(SDPX),std(SDPX),'linestyle','none','color','k');
% boxplot(SDPX,'colors','k','symbol','o')
set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('Power variability'),ylabel('(W)'); %xtickangle(45
xlim([0.5 2.5]);

%% Scatterplots correlate metrics above to performance %%%%%%%%%%%%%%%%%%%%%
% IP Force vs Sway Reduction
subplot(2,6,7),hold on
for i = 1:length(subj_array_force) % must plot one at a time to get diff color point
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanFx(i,2),SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
% First check if data is normal
% [r(1),p(1)] = corr(meanFx(:,2),SwayVRedF)
[p,r,test] = corr2vars(meanFx(:,2),SwayVRedF)
%p = 6*p; % Bonferroni correction
if p < 0.05
    c = polyfit(meanFx(:,2),SwayVRedF,1);
    plot(meanFx(:,2),polyval(c,meanFx(:,2)),'k--');
    titlename = sprintf('Sway Red. vs. Mean Force \np = %.2f, rho = %.2f',p,r); 
else
    titlename = sprintf('Sway Red. vs. Mean Force \np = %.2f',p(1)); 
end
title(titlename),box off;
% ylim([0 0.08]),xlim([0 8]);
xlabel('(N)');ylabel('SD ML torso disp. (m)');

% IP Force Var
subplot(2,6,8),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(SDFx(i,2),SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
[p,r,test] = corr2vars(SDFx(:,2),SwayVRedF)
%p = 6*p; % Bonferroni correction
if p < 0.05
    c = polyfit(SDFx(:,2),SwayVRedF,1);
    plot(SDFx(:,2),polyval(c,SDFx(:,2)),'k--');
    titlename = sprintf('Sway Red. vs. Force Variab. \np = %.2f, rho = %.2f',p,r);
else
    titlename = sprintf('Sway Red. vs. Force Variab. \np = %.2f',p); 
end
title(titlename),box off; 
% ylim([0 0.08]),xlim([0 0.6]);
xlabel('(N)');

% Mean Vel IP 
subplot(2,6,9),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanVx(i,2),SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
[p,r,test] = corr2vars(meanVx(:,2),SwayVRedF)
%p = 6*p; % Bonferroni correction
if p < 0.05
    c = polyfit(meanVx(:,2),SwayVRedF,1);
    plot(meanVx(:,2),polyval(c,meanVx(:,2)),'k--');
    titlename = sprintf('Sway Red. vs. Mean Vel \np = %.2f, rho = %.2f',p,r);
else
    titlename = sprintf('Sway Red. vs. Mean Vel \np = %.2f',p); 
end
title(titlename),box off; 
% ylim([0 0.08]),xlim([0 0.6]);
xlabel('(m/s)');

% Vel IP Var 
subplot(2,6,10),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(SDVx(i,2),SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
[p,r,test] = corr2vars(SDVx(:,2),SwayVRedF)
%p = 6*p; % Bonferroni correction
if p < 0.05
    c = polyfit(SDVx(:,2),SwayVRedF,1);
    plot(SDVx(:,2),polyval(c,SDVx(:,2)),'k--');
    titlename = sprintf('Sway Red. vs. Vel Variab. \np = %.2f, rho = %.2f',p,r);
else
    titlename = sprintf('Sway Red. vs. Vel Variab. \np = %.2f',p); 
end
title(titlename),box off; 
% ylim([0 0.08]),xlim([0 0.6]);
xlabel('(m/s)');

% Mean P IP 
subplot(2,6,11),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanPX(i,2),SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
[p,r,test] = corr2vars(meanPX(:,2),SwayVRedF)
%p = 6*p; % Bonferroni correction
if p < 0.05
    c = polyfit(meanPX(:,2),SwayVRedF,1);
    plot(meanPX(:,2),polyval(c,meanPX(:,2)),'k--');
    titlename = sprintf('Sway Red. vs. Mean Power \np = %.2f, rho = %.2f',p,r);
else
    titlename = sprintf('Sway Red. vs. Mean Power \np = %.2f',p); 
end
title(titlename),box off; 
% ylim([0 0.08]),xlim([0 0.6]);
xlabel('(W)');

% P IP Var 
subplot(2,6,12),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(SDPX(i,2),SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
[p,r,test] = corr2vars(SDPX(:,2),SwayVRedF)
%p = 6*p; % Bonferroni correction
if p < 0.05
    c = polyfit(SDPX(:,2),SwayVRedF,1);
    plot(SDPX(:,2),polyval(c,SDPX(:,2)),'k--');
    titlename = sprintf('Sway Red. vs. Power Variab. \np = %.2f, rho = %.2f',p,r);
else
    titlename = sprintf('Sway Red. vs. Power Variab. \np = %.2f',p); 
end
title(titlename),box off; 
% ylim([0 0.08]),xlim([0 0.6]);
xlabel('(W)');

%% Jitter noise correlation plots using force
% Plot each participant with error bars of +/- 5.7N

% subplot(2,6,7),
hold on
for i = 1:length(subj_array_force) % must plot one at a time to get diff color point
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanFx(i,2),SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
    errorbar(meanFx(i,2),SwayVRedF(i),5.7*ones(size(meanFx(i,2))),'horizontal','markersize',14,'color',colors(indC,:))
end
[r,p] = corr2vars(meanFx(:,2),SwayVRedF)
%p = 6*p; % Bonferroni correction
if p < 0.05
    c = polyfit(meanFx(:,2),SwayVRedF,1);
    plot(meanFx(:,2),polyval(c,meanFx(:,2)),'k--');
    titlename = sprintf('Sway Red. vs. Mean Force \np = %.2f, rho = %.2f',p,r); 
else
    titlename = sprintf('Sway Red. vs. Mean Force \np = %.2f',p); 
end
title(titlename),box off;
% ylim([0 0.08]),xlim([0 8]);
xlabel('(N)');ylabel('SD ML torso disp. (m)');

%% Plot regression POB color coded

for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    
    subplot(1,4,1)
    hold on;
    plot(1,meanRsqxPOB(n,2),'.','markersize',14,'color',colors(ind,:));
    
    subplot(1,4,2)
    hold on;
    plot(1,meanMxPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,3)
    hold on;
    plot(1,meanBxPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,4)
    hold on;
    plot(1,meanKxPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
end

plotind = 0;

% ML dir
plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanRsqxPOB(:,2),'colors','k','symbol','o'); hold on;
title('R^2'); box off; set(gca,'xtick',5);

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanMxPOB(:,2),'colors','k','symbol','o'); 
title('Mass'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(kg)')

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanBxPOB(:,2),'colors','k','symbol','o');
title('Damping'); hline(0,'k--');set(gca,'xtick',5); box off; 
set(gca,'tickdir','out');
ylabel('(Ns/m)');

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanKxPOB(:,2),'colors','k','symbol','o');
title('Stiffness'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(N/m)')

s = [112   602   854   258];
set(gcf,'outerposition',s)

%% Plot regression POB standardized coeff's, color coded

for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    
    subplot(1,4,1)
    hold on;
    plot(1,meanRsqxPOBst(n,2),'.','markersize',14,'color',colors(ind,:));
    
    subplot(1,4,2)
    hold on;
    plot(1,meanMxPOBst(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,3)
    hold on;
    plot(1,meanBxPOBst(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,4)
    hold on;
    plot(1,meanKxPOBst(n,2),'.','markersize',14,'color',colors(ind,:)); 
end

plotind = 0;

% ML dir
plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanRsqxPOBst(:,2),'colors','k','symbol','o'); hold on;
title('R^2'); box off; set(gca,'xtick',5);

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanMxPOBst(:,2),'colors','k','symbol','o'); 
title('Mass'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(kg)')

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanBxPOBst(:,2),'colors','k','symbol','o');
title('Damping'); hline(0,'k--');set(gca,'xtick',5); box off; 
set(gca,'tickdir','out');
ylabel('(Ns/m)');

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanKxPOBst(:,2),'colors','k','symbol','o');
title('Stiffness'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(N/m)')

s = [112   602   854   258];
set(gcf,'outerposition',s)

%% Correlations of model to metrics
numrows = 2; numcols = 4;
% Correlation of model to sway reduction in x dir

% mass
ind = 1;
subplot(numrows,numcols,ind),hold on
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanMxPOB(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Mx','(N/(m/s^2))','SD torso_x disp. (m)'); 

% Damping (b)
ind = ind + 1;
subplot(numrows,numcols,ind),hold on
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanBxPOB(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Bx','(N/(m/s))',[]); 

% Stiffness (k)
ind = ind + 1;
subplot(numrows,numcols,ind),hold on
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanKxPOB(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Kx','(N/m)',[]); 

% R^2
ind = ind + 1;
subplot(numrows,numcols,ind),hold on
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanRsqxPOB(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. R^2',[],[]); 

% Correlations of model to mean Fx

% Mass (m)
ind = ind + 1;
subplot(numrows,numcols,ind),hold on
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanRsqxPOB(:,2),meanFx(:,2),subj_array,subj_array_force,colors,'Mean Fx vs. Mx','(N/(m/s^2))','SD torso_x disp. (m)'); 

% Damping (b)
ind = ind + 1;
subplot(numrows,numcols,ind),hold on
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanBxPOB(:,2),meanFx(:,2),subj_array,subj_array_force,colors,'Mean Fx vs. Bx','(N/(m/s))',[]); 

% Stiffness (k)
ind = ind + 1;
subplot(numrows,numcols,ind),hold on
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanKxPOB(:,2),meanFx(:,2),subj_array,subj_array_force,colors,'Mean Fx vs. Kx','(N/m)',[]); 

% R^2
ind = ind + 1;
subplot(numrows,numcols,ind),hold on
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanRsqxPOB(:,2),meanFx(:,2),subj_array,subj_array_force,colors,'Mean Fx vs. R^2',[],[]); 

%% Fz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot means and SD's for IP force, vel, and power Overground vs. Beam z dir
% bar plots
% figure
condLab{1} = 'Ground';
condLab{2} = 'Beam';

% Means and SD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean Force IP
subplot(2,6,1),hold on;
h = barPlot2cond(meanFz,condLab,'Mean Fz','Tension > 0, (N)');

% SD Force IP
subplot(2,6,2),hold on;
h = barPlot2cond(SDFz,condLab,'Fz variability','(N)');
sigstar({[1,2]});

% Mean Vel IP
subplot(2,6,3),hold on;
h = barPlot2cond(meanVz,condLab,'Mean RMS Vz','(m/s)');

% SD Vel IP
subplot(2,6,4),hold on;
h = barPlot2cond(SDVz,condLab,'SD RMS Vz','(m/s)');

% Mean Power IP
subplot(2,6,5),hold on;
h = barPlot2cond(meanPZ,condLab,'Mean Pz','(W)');
sigstar({[1,2]});

% SD Power IP
subplot(2,6,6),hold on;
h = barPlot2cond(SDPZ,condLab,'Pz Variability','(W)');

% Scatterplots correlate metrics above to sway reduction %%%%%%%%%%%%%%%%%%
% Calculate stats and plot here
clear p r test ind

% IP Force 
subplot(2,6,7),hold on
ind = 1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanFz(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Mean Fz','(N)','SD torso_x disp. (m)'); 

% IP Force Var
subplot(2,6,8),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(SDFz(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Fz Variab.','(N)',[]);

% Mean Vel IP 
subplot(2,6,9),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanVz(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Mean Vz','(m/s)',[]);

% Vel IP Var 
subplot(2,6,10),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(SDVx(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Vz Variab.','(m/s)',[]);

% Mean P IP 
subplot(2,6,11),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanPZ(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. mean Pz','(W)',[]);

% P IP Var 
subplot(2,6,12),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(SDPZ(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Pz Variab.','(W)',[]);

%% Plot regression POB color coded

for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    
    subplot(1,4,1)
    hold on;
    plot(1,meanRsqzPOB(n,2),'.','markersize',14,'color',colors(ind,:));
    
    subplot(1,4,2)
    hold on;
    plot(1,meanMzPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,3)
    hold on;
    plot(1,meanBzPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,4)
    hold on;
    plot(1,meanKzPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
end

plotind = 0;

% ML dir
plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanRsqzPOB(:,2),'colors','k','symbol','o'); hold on;
title('R^2'); box off; set(gca,'xtick',5);

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanMzPOB(:,2),'colors','k','symbol','o'); 
title('Mass'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('N/(m/s^2)')

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanBzPOB(:,2),'colors','k','symbol','o');
title('Damping'); hline(0,'k--');set(gca,'xtick',5); box off; 
set(gca,'tickdir','out');
ylabel('N/(m/s)');

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanKzPOB(:,2),'colors','k','symbol','o');
title('Stiffness'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('N/m')

s = [112   602   854   258];
set(gcf,'outerposition',s)

%% Fxz 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot means and SD's for IP force, vel, and power Overground vs. Beam xz plane
% bar plots
% figure
condLab{1} = 'Ground';
condLab{2} = 'Beam';

% Means and SD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean Force IP
subplot(2,6,1),hold on;
h = barPlot2cond(meanFxz,condLab,'Mean Fxz','Tension > 0, (N)');
sigstar({[1,2]});

% SD Force IP
subplot(2,6,2),hold on;
h = barPlot2cond(SDFxz,condLab,'Fxz variability','(N)');
sigstar({[1,2]});

% Mean Theta 
subplot(2,6,3),hold on;
h = barPlot2cond(meanTheta,condLab,'Mean theta','deg');

% SD Theta
subplot(2,6,4),hold on;
h = barPlot2cond(SDtheta,condLab,'SD theta','deg');
sigstar({[1,2]});

% % Mean Vel IP
% subplot(2,6,3),hold on;
% h = barPlot2cond(meanVz,condLab,'Mean RMS Vz','(m/s)');
% 
% % SD Vel IP
% subplot(2,6,4),hold on;
% h = barPlot2cond(SDVz,condLab,'SD RMS Vz','(m/s)');

% Mean Power IP
subplot(2,6,5),hold on;
h = barPlot2cond(meanPXZ,condLab,'Mean Pxz','(W)');

% SD Power IP
subplot(2,6,6),hold on;
h = barPlot2cond(SDPXZ,condLab,'Pxz Variability','(W)');

% Scatterplots correlate metrics above to sway reduction %%%%%%%%%%%%%%%%%%
% Calculate stats and plot here
clear p r test ind

% IP Force 
subplot(2,6,7),hold on
ind = 1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanFxz(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Mean Fxz','(N)','SD torso_z disp. (m)'); 

% IP Force Var
subplot(2,6,8),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(SDFxz(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Fxz Variab.','(N)',[]);

% Mean Vel IP 
subplot(2,6,9),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanTheta(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Mean Theta','(deg)',[]);

% Vel IP Var 
subplot(2,6,10),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(SDtheta(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Theta Variab.','(deg)',[]);

% Mean P IP 
subplot(2,6,11),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanPXZ(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. mean Pxz','(W)',[]);

% P IP Var 
subplot(2,6,12),hold on
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(SDPXZ(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Pxz Variab.','(W)',[]);

%% Plot force vectors along theta of magnitude Fxz for all participants during beam-walking

hold on;
for cond = 1:2
    subplot(1,2,cond),hold on;
    for n = 1:length(meanTheta)
        subj = subj_array_force(n);
        thetaRad = meanTheta(n,cond)*pi/180;
        u = meanFxz(n,cond)*cos(thetaRad);
        v = meanFxz(n,cond)*sin(thetaRad);

        indC = find(subj_array == subj,1,'first'); % Use same color as kinem data

        h = quiver(0,0,u,v,'color',colors(indC,:));
    end
    axis square
    xlabel('Fxz*cos(theta)'),ylabel('Fxz*sin(theta)')
    if cond == 1
        title('Overground');
    else
        title('Beam');
    end
end

%% Torque %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bar plot means for TyGd, wTorso, PangGd Overground vs. Beam
% Means only (SD of RMS value is weird)
figure

numrows = 2; numcols = 5;
condLab{1} = 'Ground';
condLab{2} = 'Beam';

plotind = 0;
% Mean Ty about ground axis
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanTyGd,condLab,'Mean RMS Ty gd','(Nm)');
sigstar({[1,2]});

% wTorso
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanTorsoW,condLab,'RMS Torso ang speed (rad/s)','');
sigstar({[1,2]});

% Angular power RMS
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanPangGdRMS,condLab,'RMS Angular Power (N*m/s)','');

% Pos angular power 
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanPangGdPos,condLab,'Pos Ang Power (N*m/s)','');

% Neg angular power 
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanPangGdNeg,condLab,'Neg Ang Power (N*m/s)','');

% % Angular power signed
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind),hold on;
% h = barPlot2cond(meanPangGd,condLab,'Angular Power (N*m/s)','');
% sigstar({[1,2]});
% 

% Scatterplots correlate metrics above to sway reduction %%%%%%%%%%%%%%%%%%
% Calculate stats and plot here
clear p r test ind

% Ty 
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
ind = 1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanTyGd(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Mean Ty gd','(Nm)','SD ML torso disp. (m)'); 

% wTorso
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanTorsoW(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. RMS Torso ang speed','(rad/s)',[]);

% Angular power RMS
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanPangGdRMS(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. RMS Angular Power','(N*m/s)',[]);

% Pos angular power
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanPangGdPos(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Pos Ang Power','(N*m/s)',[]);

% Neg angular power
plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
ind = ind+1;
[h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanPangGdNeg(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Neg Ang Power','(N*m/s)',[]);

% % Angular power
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind),hold on;
% ind = ind+1;
% [h, p(ind), r(ind), test{ind}] = plotBivarCorr(meanPangGd(:,2),SwayVRedF,subj_array,subj_array_force,colors,'Sway Red. vs. Angular Power','(N*m/s)',[]);

%% Plot xcorr Ty to torso state, informs lag for regression
% Torque lags torso state variable if lag < 0

figure
numrows = 2; numcols = 3;
condLab{1} = 'Ground';
condLab{2} = 'Beam';

plotind = 0;

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanLagTyAngTorso,condLab,'Mean lag Ty angTorso','(s)');

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanLagTyWTorso,condLab,'Mean lag Ty wTorso','(s)');

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanLagTyAlphaTorso,condLab,'Mean lag Ty alphaTorso','(s)');

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanXcorrTyAngTorso,condLab,'Mean xcorr Ty angTorso','');

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanXcorrTyWTorso,condLab,'Mean xcorr Ty wTorso','');

plotind = plotind + 1;
subplot(numrows,numcols,plotind),hold on;
h = barPlot2cond(meanXcorrTyAlphaTorso,condLab,'Mean xcorr Ty alphaTorso','');

%% Plot regression POB angular color coded (Partnered beam cond)
% Also calculate thin rod inertia and plot to compare vs. fit inertia
numrows = 1; numcols = 4;
ind = 0;
for subj = subj_array_force
    ind = ind + 1;
    % Get mass of subj for ang momentum metric. Since compare solo to
    % assist beam, only care subj's with dynamics data
    if subj == 3
        mass = 49.3; % kg
        height = 1.60; % m
    elseif subj == 4
        mass = 56.9; % kg
        height = 1.61; % m
    elseif subj == 5
        mass = 73.8; % kg
        height = 1.64; % m
    elseif subj == 6
        mass = 78.2; % kg
        height = 1.73; % m
    elseif subj == 7
        mass = 76.0; % kg
        height = 1.82; % m
    elseif subj == 8
        mass = 67.5; % kg
        height = 1.69; % m
    elseif subj == 9
        mass = 75.2; % kg
        height = 1.68; % m
    elseif subj == 10
        mass = 55.8; % kg
        height = 1.61; % m
    elseif subj == 11
        mass = 62.0; % kg
        height = 1.64; % m
    elseif subj == 12
        mass = 67.2; % kg
        height = 1.72; % m
    elseif subj == 13
        mass = 54.8; % kg
        height = 1.61; % m
    elseif subj == 14
        mass = 69.3; % kg
        height = 1.78; % m
    end
    Icalc(ind) = mass*height^2/3;
end 

for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    
    subplot(numrows,numcols,1)
    hold on;
    plot(1,meanRsqangPOB(n,2),'.','markersize',14,'color',colors(ind,:));
    
    subplot(numrows,numcols,2)
    hold on;
    plot(1,meanMangPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
%     subplot(numrows,numcols,3)
%     hold on;
%     plot(1,Icalc(n),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(numrows,numcols,3)
    hold on;
    plot(1,meanBangPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(numrows,numcols,4)
    hold on;
    plot(1,meanKangPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
end

plotind = 0;

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(meanRsqangPOB(:,2),'colors','k','symbol','o'); hold on;
title('R^2'); box off; set(gca,'xtick',5);

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(meanMangPOB(:,2),'colors','k','symbol','o'); 
title('Inertia fit'); hline(0,'k--');
box off; set(gca,'xtick',5); % set(gca,'xtick',1:2,'xticklabel',{'Fit','Est'});
ylabel('(kg*m^2)')

% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(Icalc,'colors','k','symbol','o');
% title('Inertia est'); hline(0,'k--');set(gca,'xtick',5); box off; 
% set(gca,'tickdir','out');
% ylabel('(kg*m^2)');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(meanBangPOB(:,2),'colors','k','symbol','o');
title('Damping'); hline(0,'k--');set(gca,'xtick',5); box off; 
set(gca,'tickdir','out');
ylabel('(N*m)/(rad/s)');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(meanKangPOB(:,2),'colors','k','symbol','o');
title('Stiffness'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(N*m)/rad)')

s = [112   602   854   258];
set(gcf,'outerposition',s)

%% Plot regression POB angular color coded (Partnered ground cond)
% Also calculate thin rod inertia and plot to compare vs. fit inertia
numrows = 1; numcols = 4;
ind = 0;
for subj = subj_array_force
    ind = ind + 1;
    % Get mass of subj for ang momentum metric. Since compare solo to
    % assist beam, only care subj's with dynamics data
    if subj == 3
        mass = 49.3; % kg
        height = 1.60; % m
    elseif subj == 4
        mass = 56.9; % kg
        height = 1.61; % m
    elseif subj == 5
        mass = 73.8; % kg
        height = 1.64; % m
    elseif subj == 6
        mass = 78.2; % kg
        height = 1.73; % m
    elseif subj == 7
        mass = 76.0; % kg
        height = 1.82; % m
    elseif subj == 8
        mass = 67.5; % kg
        height = 1.69; % m
    elseif subj == 9
        mass = 75.2; % kg
        height = 1.68; % m
    elseif subj == 10
        mass = 55.8; % kg
        height = 1.61; % m
    elseif subj == 11
        mass = 62.0; % kg
        height = 1.64; % m
    elseif subj == 12
        mass = 67.2; % kg
        height = 1.72; % m
    elseif subj == 13
        mass = 54.8; % kg
        height = 1.61; % m
    elseif subj == 14
        mass = 69.3; % kg
        height = 1.78; % m
    end
    Icalc(ind) = mass*height^2/3;
end 

figure
for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    
    subplot(numrows,numcols,1)
    hold on;
    plot(1,meanRsqangPOB(n,1),'.','markersize',14,'color',colors(ind,:));
    
    subplot(numrows,numcols,2)
    hold on;
    plot(1,meanMangPOB(n,1),'.','markersize',14,'color',colors(ind,:)); 
    
%     subplot(numrows,numcols,3)
%     hold on;
%     plot(1,Icalc(n),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(numrows,numcols,3)
    hold on;
    plot(1,meanBangPOB(n,1),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(numrows,numcols,4)
    hold on;
    plot(1,meanKangPOB(n,1),'.','markersize',14,'color',colors(ind,:)); 
end

plotind = 0;

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(meanRsqangPOB(:,1),'colors','k','symbol','o'); hold on;
title('R^2'); box off; set(gca,'xtick',5);

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(meanMangPOB(:,1),'colors','k','symbol','o'); 
title('Inertia fit'); hline(0,'k--');
box off; set(gca,'xtick',5); % set(gca,'xtick',1:2,'xticklabel',{'Fit','Est'});
ylabel('(kg*m^2)')

% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(Icalc,'colors','k','symbol','o');
% title('Inertia est'); hline(0,'k--');set(gca,'xtick',5); box off; 
% set(gca,'tickdir','out');
% ylabel('(kg*m^2)');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(meanBangPOB(:,1),'colors','k','symbol','o');
title('Damping'); hline(0,'k--');set(gca,'xtick',5); box off; 
set(gca,'tickdir','out');
ylabel('(N*m)/(rad/s)');

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
boxplot(meanKangPOB(:,1),'colors','k','symbol','o');
title('Stiffness'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(N*m)/rad)')

s = [112   602   854   258];
set(gcf,'outerposition',s)

%% Old plots after this point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlations of sway (not improvement in balance) to Fx

% IP Force vs Sway
[c,ia,ib] = intersect(subj_array,subj_array_force);
SwayF = StdSway(ia,:);
figure
hold on;
for i = 1:length(subj_array_force) % must plot one at a time to get diff color point
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(SwayF(i,4),meanFx(i,2),'.','markersize',14,'color',colors(indC,:));
end
[r,p] = corr(SwayF(:,4),meanFx(:,2))

if p < 0.05
    c = polyfit(SwayF(i,4),meanFx(i,2),1);
    plot(meanFx(:,2),polyval(c,SwayF(i,4)),'k--');
    titlename = sprintf('Sway vs. Mean Force \np = %.2f, rho = %.2f',p,r); 
else
    titlename = sprintf('Sway vs. Mean Force \np = %.2f',p); 
end
title(titlename),box off;
xlabel('(m)');ylabel('(N)');
%% % Power POB torso in x dir
% subplot(2,3,3),hold on;
% for n = 1:length(subj_array_force)
%     subj = subj_array_force(n);
%     % Find color
%     ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(meanPmagPOBX(n,:),'-','color',colors(ind,:));
% end
% boxplot(meanPmagPOBX,'colors','k','symbol','o'); sigstar({[1,2]}); 
% xlim([0.5 2.5]);
% set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
% title('Mean power POB magnitude'),ylabel('(W)'); %xtickangle(45

% % SD's %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Force
% subplot(2,3,4),hold on;
% for n = 1:length(subj_array_force)
%     subj = subj_array_force(n);
%     % Find color
%     ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(SDFx(n,:),'-','color',colors(ind,:));
% end
% boxplot(SDFx,'colors','k','symbol','o'); sigstar({[1,2]}); 
% set(gca,'xticklabel',condLab); box off; set(gca,'tickdir','out');
% title('Force magnitude variability'),ylabel('(N)'); %xtickangle(45)
% xlim([0.5 2.5]);

% % Power IP
% subplot(2,3,5),hold on;
% for n = 1:length(subj_array_force)
%     subj = subj_array_force(n);
%     % Find color
%     ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(SDPmagX(n,:),'-','color',colors(ind,:));
% end
% boxplot(SDPmagX,'colors','k','symbol','o');
% xlim([0.5 2.5]);
% set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
% title('Power IP magnitude variability'),ylabel('(W)');
% 
% % Power POB torso
% subplot(2,3,6),hold on;
% for n = 1:length(subj_array_force)
%     subj = subj_array_force(n);
%     % Find color
%     ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(SDPmagPOBX(n,:),'-','color',colors(ind,:));
% end
% boxplot(SDPmagPOBX,'colors','k','symbol','o');
% xlim([0.5 2.5]);
% set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
% title('Power POB magnitude variability'),ylabel('(W)');

%% Compare positive vs. negative force magnitude for Assist Beam
% Some participants don't have any negative forces for some trials, so just
% plot a '.' for their positive value

% Plot color coded

% Force
xlab{1} = 'Push';
xlab{2} = 'Pull';
subplot(1,3,1),hold on;
% Plot each subj with color
for n = 1:length(subj_array_force)
    subj = subj_array_force(n)
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    if isnan(meanNegFx(n,2)) % no neg force
        plot(2,meanPosFx(n,2),'.','markersize',14,'color',colors(ind,:));
    else
        plot([-meanNegFx(n,2) meanPosFx(n,2)],'-','color',colors(ind,:));
    end
end
% Box plot
boxplot([-meanNegFx(:,2) meanPosFx(:,2)],'colors','k','symbol','o');
sigstar({[1,2]});
xlim([0.5 2.5]);
set(gca,'xticklabel',xlab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('Force magnitude'),ylabel('(N)');

% Force - worst-case scenario for drift
meanPosFx(:,2) = meanPosFx(:,2) - 1.5;
meanNegFx(:,2) = meanNegFx(:,2) - 1.5;
xlab{1} = 'Push';
xlab{2} = 'Pull';
subplot(1,3,2),hold on;
% Plot each subj with color
for n = 1:length(subj_array_force)
    subj = subj_array_force(n)
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    if isnan(meanNegFx(n,2)) % no neg force
        plot(2,meanPosFx(n,2),'.','markersize',14,'color',colors(ind,:));
    else
        plot([-meanNegFx(n,2) meanPosFx(n,2)],'-','color',colors(ind,:));
    end
end
% Box plot
boxplot([-meanNegFx(:,2) meanPosFx(:,2)],'colors','k','symbol','o');
sigstar({[1,2]});
xlim([0.5 2.5]);
set(gca,'xticklabel',xlab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('Force magnitude'),ylabel('(N)');

%% Power IP pos vs. neg
xlab{1} = 'Motor';
xlab{2} = 'Brake';
% subplot(1,3,1),
hold on
for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot([-meanPnegX(n,2) meanPposX(n,2)],'-','color',colors(ind,:)); % pos power is brake
end
% Box plot
boxplot([-meanPnegX(:,2) meanPposX(:,2)],'colors','k','symbol','o'); 
sigstar({[1,2]});
xlim([0.5 2.5]);
set(gca,'xticklabel',xlab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('IP power'),ylabel('(W)');

%% Power pos and neg from signProdPer
figure
xlab{1} = 'Motor';
xlab{2} = 'Brake';
subplot(1,3,1),hold on
for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot([-meanPnegF(n,2) meanPposF(n,2)],'-','color',colors(ind,:)); % pos power is brake
end
% Box plot
boxplot([-meanPnegF(:,2) meanPposF(:,2)],'colors','k','symbol','o'); 
sigstar({[1,2]});
xlim([0.5 2.5]);
set(gca,'xticklabel',xlab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('IP power'),ylabel('(W)');

% Worst case F drift lo
subplot(1,3,2),hold on
for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot([-meanPnegFlo(n,2) meanPposFlo(n,2)],'-','color',colors(ind,:)); % pos power is brake
end
% Box plot
boxplot([-meanPnegFlo(:,2) meanPposFlo(:,2)],'colors','k','symbol','o'); 
sigstar({[1,2]});
xlim([0.5 2.5]);
set(gca,'xticklabel',xlab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('IP power F drift low bound'),ylabel('(W)');

% Worst case F drift hi
subplot(1,3,3),hold on
for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot([-meanPnegFhi(n,2) meanPposFhi(n,2)],'-','color',colors(ind,:)); % pos power is brake
end
% Box plot
boxplot([-meanPnegFhi(:,2) meanPposFhi(:,2)],'colors','k','symbol','o'); 
sigstar({[1,2]});
xlim([0.5 2.5]);
set(gca,'xticklabel',xlab,'xtick',1:2); box off; set(gca,'tickdir','out');
title('IP power F drift hi bound'),ylabel('(W)');

% % Power POB
% xlab{1} = 'In dir';
% xlab{2} = 'Opp dir';
% subplot(1,3,3),hold on
% for n = 1:length(subj_array_force)
%     subj = subj_array_force(n);
%     % Find color
%     ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot([-meanPnegPOBX(n,2) meanPposPOBX(n,2)],'-','color',colors(ind,:)); % pos power is brake
% end
% % Box plot
% boxplot([-meanPnegPOBX(:,2) meanPposPOBX(:,2)],'colors','k','symbol','o'); 
% sigstar({[1,2]});
% xlim([0.5 2.5]);
% set(gca,'xticklabel',xlab,'xtick',1:2); box off; set(gca,'tickdir','out');
% title('Power POB magnitude'),ylabel('(W)');

%% Plot each trial's percent time spent in quadrant of F and P_IP combo
hold on;
% use same colors as kinem data with 12 subj's
plotOpt = 3; % 1: plot distances on axes, 2: plot size of circle in middle 
% of quadrant scaled to proportion of trials present, 3: plot one point per
% trial, 4: plot mean over all trials at mean with size of cirlce proportional
% to SD over all trials for each quadrant separately, 5: similar to 3 but
% plot 4 pts per trial
n = 0;
for subj = subj_array_force
    perPF = [];
    n = n + 1;
    legnames{n} = sprintf('Partnership %i',n);
    rows = GroupStats.Subject==subj & strcmp(GroupStats.Type,'Assist Beam');
    perPF{1} = GroupStats.perPposFposX(rows); % Q1
    perPF{4} = GroupStats.perPposFnegX(rows); % Q4
    perPF{2} = GroupStats.perPnegFposX(rows); % Q2
    perPF{3} = GroupStats.perPnegFnegX(rows);
    perPpos = GroupStats.perPposFposX(rows) + GroupStats.perPposFnegX(rows);
    perFpos = GroupStats.perPposFposX(rows) + GroupStats.perPnegFposX(rows);
    perPneg = GroupStats.perPnegFposX(rows) + GroupStats.perPnegFnegX(rows);
    perFneg = GroupStats.perPposFnegX(rows) + GroupStats.perPnegFnegX(rows);
    % Replace zeros with nan's bc don't want to plot those
    for i = 1:4
        a = perPF{i};
        a(a==0) = nan;
        perPF{i} = a;
    end

    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    
    if plotOpt == 1
        plot(perPF{1},perPF{1},'.','color',colors(ind,:)); % QI
        plot(perPF{4},-perPF{4},'.','color',colors(ind,:)); % Q1V
        plot(-perPF{2},perPF{2},'.','color',colors(ind,:)); % QII
        plot(-perPF{3},-perPF{3},'.','color',colors(ind,:));
    elseif plotOpt == 2
        posQ(1,:) = [0.5 0.5];
        posQ(2,:) = [-0.5 0.5];
        posQ(3,:) = [-0.5 -0.5];
        posQ(4,:) = [0.5 -0.5];
        for i = 1:4 % all 4 quadrants
            d(i) = nansum(perPF{i})/length(perPF{1}); % normalize sum of all % time spent data so max diameter of circle is 1
            r = rectangle('Position',[posQ(i,:)-d(i)/2*[1 1] d(i) d(i)],'Curvature',1); % plot circle centered about posQ1
            set(r,'facecolor',[colors(ind,:) 0.5]); % 4th element of color gives opacity
        end
    elseif plotOpt == 3
        plot(perPpos,perFpos,'.','markersize',14,'color',[colors(ind,:) 0.1]);
    elseif plotOpt == 4
        plotCircle(perPFmean(n),perPFmean(n),perPFSD(n));
    elseif plotOpt == 5
        plot(perPpos,perFpos,'.','markersize',14,'color',[colors(ind,:) 0.1]); % Q1
        plot(-perPneg,perFpos,'.','markersize',14,'color',[colors(ind,:) 0.1]); % Q2
        plot(-perPneg,-perFneg,'.','markersize',14,'color',[colors(ind,:) 0.1]); % Q3
        plot(perPpos,-perFneg,'.','markersize',14,'color',[colors(ind,:) 0.1]); % Q4
    end
end
if plotOpt == 3
    legend(legnames);
end
% xlim([-1 1]),ylim([-1 1]); 
xlim([0 1]),ylim([0 1]); 
axis square
set(gca,'xtick',2:4,'ytick',2:4) % no ticks
% vline(0,'k-');hline(0,'k-'); box on;
vline(0.5,'k-');hline(0.5,'k-'); box on;

%% Plot each trial's percent time spent in quadrant of oppose/amplify deviation/return
 
hold on;
% use same colors as kinem data with 12 subj's
% plot one point per trial
n = 0;
for subj = subj_array_force
    perPF = [];
    n = n + 1;
    legnames{n} = sprintf('Partnership %i',n);
%     rows = GroupStats.Subject==subj & strcmp(GroupStats.Type,'Assist Beam');
    rows = GroupStats.Subject==subj & strcmp(GroupStats.Type,'Assist Ground');
    perOpp = GroupStats.POBperOppDevX(rows) + GroupStats.POBperOppRetX(rows);
    perRet = GroupStats.POBperOppRetX(rows) + GroupStats.POBperAmpRetX(rows);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data  
    plot(perOpp,perRet,'.','markersize',14,'color',[colors(ind,:) 0.1]);
end
% legend(legnames);
 
xlim([0 1]),ylim([0 1]); 
axis square
set(gca,'xtick',2:4,'ytick',2:4) % no ticks
% vline(0,'k-');hline(0,'k-'); box on;
vline(0.5,'k-');hline(0.5,'k-'); box on;

%% Plot each trial's percent time spent in quadrant of F and P_POB combo

hold on;
% use same colors as kinem data with 12 subj's
% plot one point per trial
n = 0;
for subj = subj_array_force
    n = n + 1;
    legnames{n} = sprintf('Partnership %i',n);
    rows = GroupStats.Subject==subj & strcmp(GroupStats.Type,'Assist Beam');
    perPposPOB = GroupStats.POBperPposFposX(rows) + GroupStats.POBperPposRetX(rows);
    perFpos = GroupStats.POBperPposFposX(rows) + GroupStats.POBperPnegFposX(rows);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data  
    plot(perPposPOB,perFpos,'.','markersize',14,'color',[colors(ind,:) 0.1]);
end
% legend(legnames);
 
xlim([0 1]),ylim([0 1]); 
axis square
set(gca,'xtick',2:4,'ytick',2:4) % no ticks
% vline(0,'k-');hline(0,'k-'); box on;
vline(0.5,'k-');hline(0.5,'k-'); box on;

%% Plot regression IP color coded

for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    
    subplot(1,4,1)
    hold on;
    plot(1,meanRsqxIP(n,2),'.','markersize',14,'color',colors(ind,:));
    
    subplot(1,4,2)
    hold on;
    plot(1,meanMxIP(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,3)
    hold on;
    plot(1,meanBxIP(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,4)
    hold on;
    plot(1,meanKxIP(n,2),'.','markersize',14,'color',colors(ind,:)); 
end

plotind = 0;

% ML dir
plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanRsqxIP(:,2),'colors','k','symbol','o'); hold on;
title('R^2'); box off; set(gca,'xtick',5);

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanMxIP(:,2),'colors','k','symbol','o'); 
title('Mass'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(kg)')

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanBxIP(:,2),'colors','k','symbol','o');
title('Damping'); hline(0,'k--');set(gca,'xtick',5); box off; 
set(gca,'tickdir','out');
ylabel('(Ns/m)');

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanKxIP(:,2),'colors','k','symbol','o');
title('Stiffness'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(N/m)')

s = [112   602   854   258];
set(gcf,'outerposition',s)

%% Plot regression POB with lag color coded
figure
for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    
    subplot(1,4,1)
    hold on;
    plot(1,meanRsqxLagPOB(n,2),'.','markersize',14,'color',colors(ind,:));
    
    subplot(1,4,2)
    hold on;
    plot(1,meanMxLagPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,3)
    hold on;
    plot(1,meanBxLagPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
    
    subplot(1,4,4)
    hold on;
    plot(1,meanKxLagPOB(n,2),'.','markersize',14,'color',colors(ind,:)); 
end

plotind = 0;

% ML dir
plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanRsqxLagPOB(:,2),'colors','k','symbol','o'); hold on;
title('R^2'); box off; set(gca,'xtick',5);

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanMxLagPOB(:,2),'colors','k','symbol','o'); 
title('Mass'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(kg)')

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanBxLagPOB(:,2),'colors','k','symbol','o');
title('Damping'); hline(0,'k--');set(gca,'xtick',5); box off; 
set(gca,'tickdir','out');
ylabel('(Ns/m)');

plotind = plotind + 1;
subplot(1,4,plotind)
boxplot(meanKxLagPOB(:,2),'colors','k','symbol','o');
title('Stiffness'); hline(0,'k--');
box off; set(gca,'xtick',5);
ylabel('(N/m)')

s = [112   602   854   258];
set(gcf,'outerposition',s)

%% Plot regression coeff's showing variab w/in subj and across subj

color_array = colors;
for n = 1:length(subj_array_force)
    subj = subj_array_force(n);
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    % R^2
    subplot(4,length(subj_array_force)+1,n)
    boxplot(meanRsqxArray{n},'colors',color_array(ind,:),'symbol','o');
    hold on; ylim([0 0.8]); box off;
    if n == 1
        ylabel('R^2');
    end
    % M
    subplot(4,length(subj_array_force)+1,length(subj_array_force)+1+n)
    boxplot(meanMxArray{n},'colors',color_array(ind,:),'symbol','o');
    hold on; hline(0,'k--');
    ylim([-5 3]); box off;
    if n == 1
        ylabel('mass (kg)');
    end
    % B
    subplot(4,length(subj_array_force)+1,2*(length(subj_array_force)+1)+n)
    boxplot(meanBxArray{n},'colors',color_array(ind,:),'symbol','o');
    hold on; hline(0,'k--');
    ylim([-15 50]); box off;
    if n == 1
        ylabel('damping');
    end
    % K
    subplot(4,length(subj_array_force)+1,3*(length(subj_array_force)+1)+n)
    boxplot(meanKxArray{n},'colors',color_array(ind,:),'symbol','o');
    hold on; hline(0,'k--'); box off;
    ylim([-70 130]);
    if n == 1
        ylabel('stiffness');
    end
end

% ML dir
subplot(4,length(subj_array_force)+1,length(subj_array_force)+1)
boxplot(meanRsqxPOB(:,2),'colors','k','symbol','o'); hold on; ylim([0 0.8]);
% title('R^2'); 
box off; set(gca,'xtick',5);

subplot(4,length(subj_array_force)+1,2*(length(subj_array_force)+1))
boxplot(meanMxPOB(:,2),'colors','k','symbol','o'); 
% title('Mass'); 
hline(0,'k--');
box off; set(gca,'xtick',5); ylim([-5 3]);
ylabel('(kg)')

subplot(4,length(subj_array_force)+1,3*(length(subj_array_force)+1))
boxplot(meanBxPOB(:,2),'colors','k','symbol','o');
% title('Damping'); 
hline(0,'k--');set(gca,'xtick',5); box off; 
set(gca,'tickdir','out'); ylim([-15 50]);
ylabel('(Ns/m)');

subplot(4,length(subj_array_force)+1,4*(length(subj_array_force)+1))
boxplot(meanKxPOB(:,2),'colors','k','symbol','o');
% title('Stiffness'); 
hline(0,'k--');
box off; set(gca,'xtick',5); ylim([-70 130]);
ylabel('(N/m)')

% s = [112   602   854   258];
% set(gcf,'outerposition',s)

%% % AP dir
% plotind = plotind + 1;
% subplot(3,4,plotind)
% boxplot(meanRsqyPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Rsq'); 
% title('Fit POB torso AP direction');
% 
% plotind = plotind + 1;
% subplot(3,4,plotind)
% boxplot(meanMyPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Mass (kg)'); hline(0,'k');
% 
% plotind = plotind + 1;
% subplot(3,4,plotind)
% boxplot(meanByPOB(:,2)); set(gca,'xticklabel','All partnerships');box off; set(gca,'tickdir','out');
% ylabel('Damping (Ns/m)'); hline(0,'k');
% 
% plotind = plotind + 1;
% subplot(3,4,plotind)
% boxplot(meanKyPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Stiffness (N/m)'); hline(0,'k');
% 
% % Vert dir
% plotind = plotind + 1;
% subplot(3,4,plotind)
% boxplot(meanRsqzPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Rsq');
% title('Fit POB torso vert direction');
% 
% plotind = plotind + 1;
% subplot(3,4,plotind)
% boxplot(meanMzPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Mass (kg)'); hline(0,'k');
% 
% plotind = plotind + 1;
% subplot(3,4,plotind)
% boxplot(meanBzPOB(:,2)); set(gca,'xticklabel','All partnerships');box off; set(gca,'tickdir','out');
% ylabel('Damping (Ns/m)'); hline(0,'k');
% 
% plotind = plotind + 1;
% subplot(3,4,plotind)
% boxplot(meanKzPOB(:,2)); set(gca,'xticklabel','All partnerships'); box off; set(gca,'tickdir','out');
% ylabel('Stiffness (N/m)'); hline(0,'k');

%% mean percent Partner Beam trials where m, b, k sig. ML dir

plotind = 0;
numrows = 3; numcols = 1;

plotind = plotind + 1;
subplot(numrows,numcols,plotind)
bar(mxPOBsig*100),ylabel('% trials m sig')
box off; ylim([0 100]); xlabel('Partnership');
set(gca,'tickdir','out');
title('Fit POB torso');

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

%% Plot xcorr lag different IP to POB torso signals for ML dir
xtlab{1} = conds{2};
xtlab{2} = conds{4};

subplot(2,3,1)
boxplot(meanLagFIPTorsoX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');
ylabel('lag at max xcorr (s)');
title('ML force and torso pos');

subplot(2,3,2)
boxplot(meanLagFIPvTorsoX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');
title('ML force and torso vel');

subplot(2,3,3)
boxplot(meanLagvIPvTorsoX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');
title('IP and torso vels');

subplot(2,3,4)
boxplot(meanxcorrFIPTorsoX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');
ylabel('max xcorr');

subplot(2,3,5)
boxplot(meanxcorrFIPvTorsoX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');

subplot(2,3,6)
boxplot(meanxcorrvIPvTorsoX); set(gca,'xticklabel',xtlab),hold on, hline(0,'k--');
box off; set(gca,'tickdir','out');

%% Plot POB arm len to check xcorr results
xtlab{1} = conds{1};
xtlab{2} = conds{2};

subplot(1,2,1)
boxplot(armPOBx),set(gca,'xticklabel',xtlab),hold on
sigstar({[1 2]});
box off; set(gca,'tickdir','out'),title('Distance IP to POB torso (m)');
ylabel('mean (m)');

subplot(1,2,2)
boxplot(SDarmPOBx),set(gca,'xticklabel',xtlab),hold on
sigstar({[1 2]});
box off; set(gca,'tickdir','out')
ylabel('SD (m)');

%% Plot metrics for Assist Beam cond vs solo beam distance
plotind = 1;

% Mean and SD F and P mag
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i)
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),meanFx(i,2),'.','markersize',14,'color',colors(indC,:));
end
title('Assist Beam'); box off;
xlabel('Solo beam dist (m)'),ylabel('Mean F mag (N)');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),meanPmagX(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('Mean P mag (W)');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),SDFx(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('Mean F var (N)');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),SDPmagX(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('Mean P var (W)');

% Pos and neg F and P mag's
plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),meanPosFx(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('Mean Tension (N)');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),-meanNegFx(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('Mean Compression (N)');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),meanPposX(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('Mean Braking Power (W)');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),-meanPnegX(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('Mean Motor Power (W)');

% %time spent in each F and P strategy
plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),perFpos(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('% time tension');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),perFneg(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('% time compression');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),perPpos(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('% time brake');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),perPneg(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('% time motor');

% Model fits
plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),meanRsqxPOB(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('R^2');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),meanMxPOB(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('mass (kg)');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),meanBxPOB(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('damping (N/(m/s))');

plotind = plotind+1;
subplot(4,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(soloDistF(i),meanKxPOB(i,2),'.','markersize',14,'color',colors(indC,:));
end
box off;
xlabel('Solo beam dist (m)'),ylabel('stiffness (N/m)');

%% % SD F mag
% plotind = plotind+1;
% subplot(3,4,plotind),hold on
% for i = 1:length(subj_array_force)
%     subj = subj_array_force(i);
%     indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(SDFx(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
% end
% [r(plotind),p(plotind)] = corr(SDFx(:,2),-SwayVRedF);
% if p(plotind) < 0.05
%     c = polyfit(SDFx(:,2),-SwayVRedF,1);
%     plot(SDFx,polyval(c,SDFx),'k--');
%     titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind));
% else
%     titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
% end
% ylim([0 0.08]),xlim([0 1.5]); title(titlename),box off;
% xlabel('Force Mag. Var. (N)');
% 
% % SD P mag
% plotind = plotind+1;
% subplot(3,4,plotind),hold on
% for i = 1:length(subj_array_force)
%     subj = subj_array_force(i);
%     indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(SDPmagX(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
% end
% 
% [r(plotind),p(plotind)] = corr(SDPmagX(:,2),-SwayVRedF);
% if p(plotind) < 0.05
%     c = polyfit(SDPmagX(:,2),-SwayVRedF,1);
%     plot(SDPmagX,polyval(c,SDPmagX),'k--');
%     titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind));
% else
%     titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
% end
% ylim([0 0.08]),xlim([0 0.5]); title(titlename),box off;
% xlabel('Power Mag. Var. (W)');

%% Tension force mag
plotind = plotind+1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanPosFx(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end

[r(plotind),p(plotind)] = corr(meanPosFx(:,2),-SwayVRedF);
if p(plotind) < 0.05
    c = polyfit(meanPosFx(:,2),-SwayVRedF,1);
    plot(meanPosFx(:,2),polyval(c,meanPosFx(:,2)),'k--');
    titlename = sprintf('p = %.2f, p_B = %.2f, rho = %.2f',p(plotind),p(plotind)*12,r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),xlim([0 8]); title(titlename),box off;
ylabel('Improvement in Sway Var. (m)'),xlabel('Pull Force (N)');

% Compression force mag
plotind = plotind+1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(-meanNegFx(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
ind = ~isnan(meanNegFx(:,2));
[r(plotind),p(plotind)] = corr(-meanNegFx(ind,2),-SwayVRedF(ind));
if p(plotind) < 0.05
    c = polyfit(-meanNegFx(ind,2),-SwayVRedF(ind),1);
    plot(-meanNegFx(ind,2),polyval(c,-meanNegFx(ind,2)),'k--');
    titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),xlim([0 3]); title(titlename),box off;
xlabel('Push Force (N)');

% Same dir sway force mag
plotind = plotind+1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanPposPOBX(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
ind = ~isnan(meanPposPOBX(:,2))
[r(plotind),p(plotind)] = corr(meanPposPOBX(ind,2),-SwayVRedF(ind));
if p(plotind) < 0.05
    c = polyfit(meanPposPOBX(ind,2),-SwayVRedF(ind),1);
    plot(meanPposPOBX(ind,2),polyval(c,meanPposPOBX(ind,2)),'k--');
    titlename = sprintf('p = %.2f, p_B = %.2f, rho = %.2f',p(plotind),p(plotind)*12,r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),
% xlim([0 8]); 
title(titlename),box off;
ylabel('Improvement in Sway Var. (m)'),xlabel('Power in sway dir. (W)');

% Opp dir sway force mag
plotind = plotind+1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(-meanPnegPOBX(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end
ind = ~isnan(meanPnegPOBX(:,2));
[r(plotind),p(plotind)] = corr(-meanPnegPOBX(ind,2),-SwayVRedF(ind));
if p(plotind) < 0.05
    c = polyfit(-meanNegFx(ind,2),-SwayVRedF(ind),1);
    plot(-meanPnegPOBX(ind,2),polyval(c,-meanPnegPOBX(ind,2)),'k--');
    titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),
% xlim([0 3]); 
title(titlename),box off;
xlabel('Power opposite sway (W)');

% % Braking power mag
% plotind = plotind+1;
% subplot(4,4,plotind),hold on
% for i = 1:length(subj_array_force)
%     subj = subj_array_force(i);
%     indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(-SwayVRedF(i),meanPposX(i,2),'.','markersize',14,'color',colors(indC,:));
% end
% 
% [r(plotind),p(plotind)] = corr(-SwayVRedF,meanRsqxPOB(:,2));
% if p(plotind) < 0.05
%     c = polyfit(-SwayVRedF,meanRsqxPOB(i,2),1);
%     plot(-SwayVRedF,polyval(c,-SwayVRedF),'k--');
% end
% titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind)); title(titlename),box off;
% xlabel('Mean Braking Power (W)');
% 
% % Motor power mag
% plotind = plotind+1;
% subplot(4,4,plotind),hold on
% for i = 1:length(subj_array_force)
%     subj = subj_array_force(i);
%     indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(-SwayVRedF(i),-meanPnegX(i,2),'.','markersize',14,'color',colors(indC,:));
% end
% 
% [r(plotind),p(plotind)] = corr(-SwayVRedF,meanRsqxPOB(:,2));
% if p(plotind) < 0.05
%     c = polyfit(-SwayVRedF,meanRsqxPOB(i,2),1);
%     plot(-SwayVRedF,polyval(c,-SwayVRedF),'k--');
% end
% titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind)); title(titlename),box off;
% xlabel('Mean Motor Power (W)');

% Percent of time spent in tension
plotind = plotind + 1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(perFpos(i,2)*100,-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end

[r(plotind),p(plotind)] = corr(perFpos(:,2)*100,-SwayVRedF);
if p(plotind) < 0.05
    c = polyfit(perFpos(:,2)*100,-SwayVRedF,1);
    plot(perFpos(:,2)*100,polyval(c,perFpos(:,2)*100),'k--');
    titlename = sprintf('p = %.2f, p_B = %.2f, rho = %.2f',p(plotind),p(plotind)*12,r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),xlim([0 100]); title(titlename),box off;
xlabel('% Time pulling');

% Percent of time spent power in same dir as sway
plotind = plotind + 1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(perPposPOB(i,2)*100,-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end

[r(plotind),p(plotind)] = corr(perPposPOB(:,2)*100,-SwayVRedF);
if p(plotind) < 0.05
    c = polyfit(perPposPOB(:,2)*100,-SwayVRedF,1);
    plot(perPposPOB(:,2)*100,polyval(c,perPposPOB(:,2)*100),'k--');
    titlename = sprintf('p = %.2f, p_B = %.2f, rho = %.2f',p(plotind),p(plotind)*12,r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),
% xlim([0 100]); 
title(titlename),box off;
xlabel('% Time same dir. sway');

% plotind = plotind+1;
% subplot(4,4,plotind),hold on
% for i = 1:length(subj_array_force)
%     subj = subj_array_force(i);
%     indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(SwayVRedF(i),perFneg(i,2),'.','markersize',14,'color',colors(indC,:));
% end

% % Prop of time spent braking
% plotind = plotind+1;
% subplot(3,4,plotind),hold on
% for i = 1:length(subj_array_force)
%     subj = subj_array_force(i);
%     indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(perPpos(i,2)*100,-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
% end
% 
% [r(plotind),p(plotind)] = corr(perPpos(:,2),-SwayVRedF);
% if p(plotind) < 0.05
%     c = polyfit(perPpos(:,2)*100,-SwayVRedF,1);
%     plot(perPpos(:,2)*100,polyval(c,perPpos(:,2)*100),'k--');
%     titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind));
% else
%     titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
% end
% ylim([0 0.08]),xlim([0 100]); title(titlename),box off;
% xlabel('% Time Braking');

% plotind = plotind+1;
% subplot(4,4,plotind),hold on
% for i = 1:length(subj_array_force)
%     subj = subj_array_force(i);
%     indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
%     plot(SwayVRedF(i),perPneg(i,2),'.','markersize',14,'color',colors(indC,:));
% end

% Model quality of fit
plotind = 9;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanRsqxPOB(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end

[r(plotind),p(plotind)] = corr(meanRsqxPOB(:,2),-SwayVRedF);
if p(plotind) < 0.05
    c = polyfit(meanRsqxPOB(:,2),-SwayVRedF,1);
    plot(meanRsqxPOB(:,2),polyval(c,meanRsqxPOB(:,2)),'k--');
    titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),xlim([0 1]); title(titlename),box off;
xlabel('R^2');

% Model mass parameter
plotind = plotind+1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanMxPOB(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end

[r(plotind),p(plotind)] = corr(meanMxPOB(:,2),-SwayVRedF);
if p(plotind) < 0.05
    c = polyfit(meanMxPOB(:,2),-SwayVRedF,1);
    plot(meanMxPOB(:,2),polyval(c,meanMxPOB(:,2)),'k--');
    titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),xlim([-2.5 1]); title(titlename),box off;
xlabel('Mass (kg)');

% Model damping parameter
plotind = plotind+1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanBxPOB(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end

[r(plotind),p(plotind)] = corr(meanBxPOB(:,2),-SwayVRedF);
if p(plotind) < 0.05
    c = polyfit(meanBxPOB(:,2),-SwayVRedF,1);
    plot(meanBxPOB(:,2),polyval(c,meanBxPOB(:,2)),'k--');
    titlename = sprintf('p = %.2f, p_B = %.2f, rho = %.2f',p(plotind),p(plotind)*12,r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),xlim([0 25]); title(titlename),box off;
xlabel('Damping (N/(m/s))');

% Model stiffness parameter
plotind = plotind+1;
subplot(3,4,plotind),hold on
for i = 1:length(subj_array_force)
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(meanKxPOB(i,2),-SwayVRedF(i),'.','markersize',14,'color',colors(indC,:));
end

[r(plotind),p(plotind)] = corr(meanKxPOB(:,2),-SwayVRedF);
if p(plotind) < 0.05
    c = polyfit(meanKxPOB(:,2),-SwayVRedF,1);
    plot(meanKxPOB(:,2),polyval(c,meanKxPOB(:,2)),'k--');
    titlename = sprintf('p = %.2f, rho = %.2f',p(plotind),r(plotind));
else
    titlename = sprintf('p = %.2f, p_B = %.2f',p(plotind),p(plotind)*12); 
end
ylim([0 0.08]),xlim([-10 100]); title(titlename),box off;
xlabel('Stiffness (N/m)');

a = findobj(gcf,'type','axes');
set(a,'tickdir','out','ytick',0:0.04:0.08)

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
        rows = GroupStats.Subject==subj & strcmp(GroupStats.Type,conds_perf{i});
        StdSway(n,i) = nanmean(GroupStats.StdSway(rows)); % SD torso sway
        AvgSpeed(n,i) = nanmean(GroupStats.AvgSpeed(rows)); % Average speed of walking on beam
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
        rows = GroupStats.Subject==subj & strcmp(GroupStats.Type,conds{i});
        % Frontal plane segment angles
%         StdPelvicObliq(n,i) = nanmean(GroupStats.StdPelvicObliq(rows));
        StdThoraxObliq(n,i) = nanmean(GroupStats.StdThoraxObliq(rows));
        StdLegObliq(n,i) = nanmean(GroupStats.StdLegObliq(rows));
%         rPelvicThoraxObliq(n,i) = nanmean(GroupStats.rPelvicThoraxObliq(rows));
        rLegThoraxObliq(n,i) = nanmean(GroupStats.rLegThoraxObliq(rows));
        % Sway metrics
        StdSway(n,i) = nanmean(GroupStats.StdSway(rows));
        StdCOMSway(n,i) = nanmean(GroupStats.StdCOMSway(rows));
        rFverttorso(n,i) = nanmean(GroupStats.rFverttorso(rows));
        rFlatCOM(n,i) = nanmean(GroupStats.rFlatCOM(rows));
        rFvertCOM(n,i) = nanmean(GroupStats.rFvertCOM(rows));
        rFlattorso(n,i) = nanmean(GroupStats.rFlattorso(rows));
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

%% First compare if Std torso sway is very different from Std COM sway

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
boxplot(rFverttorso(:,[2 4])); set(gca,'xticklabel',assistConds); box off; set(gca,'tickdir','out');
ylabel('Corr. torso sway to vert. F'), ylim([-1 1]);

subplot(2,2,2)
boxplot(rFlatCOM(:,[2 4])); set(gca,'xticklabel',assistConds); box off; set(gca,'tickdir','out');
ylabel('Corr. COM sway to lat. F'), ylim([-1 1]);

subplot(2,2,3)
boxplot(rFvertCOM(:,[2 4])); set(gca,'xticklabel',assistConds); box off; set(gca,'tickdir','out');
ylabel('Corr. COM sway to vert. F'), ylim([-1 1]);

subplot(2,2,4)
boxplot(rFlattorso(:,[2 4])); set(gca,'xticklabel',assistConds); box off; set(gca,'tickdir','out');
ylabel('Corr. torso sway to lat. F'), ylim([-1 1]);

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
% plotind = 0;
% numcols = 3; numrows = 1;

% % Abs power vs. thresh
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(meanPmagX(:,2)); sigstar({1}); 
% hline(0.0256,'k--'), ylim([0 0.21])
% set(gca,'xticklabel','Assist Beam'); box off; set(gca,'tickdir','out');
% ylabel('Power Mag ML (W)'); 
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(meanPmagY(:,2)); sigstar({1}); 
% hline(2.621,'k--'),ylim([0 2.5])
% set(gca,'xticklabel','Assist Beam'); box off; set(gca,'tickdir','out');
% ylabel('Power Mag AP (W)'); 
% 
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot(meanPmagX(:,2)); sigstar({1}); 
% hline(0.2597,'k--'),ylim([0 0.21]);
% set(gca,'xticklabel','Assist Beam'); box off; set(gca,'tickdir','out');
% ylabel('Power Mag Vert (W)'); 

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
% title('POB torso');
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

% xlab{1} = 'Motor';
% xlab{2} = 'Brake';
% plotind = plotind + 1;
% subplot(numrows,numcols,plotind)
% boxplot([powerNegPOBX(:,2) abs(powerPosPOBX(:,2))],'colors','k','symbol','o'); 
% sigstar({[1,2]});
% set(gca,'xticklabel',xlab); box off; set(gca,'tickdir','out');
% ylabel('Power magnitude (W)'); xtickangle(45)

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
