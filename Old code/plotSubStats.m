% plot individual subj's metrics and means of metrics using concatenated tables from
% LateGroup. 
clear all; clc; close all;

subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)

%% Kinematics
kinem = load('HHI2017_LateStats_MW');
conds = {'Solo Beam','Assist Ground','Assist Beam'};
numrows = 4; numcols = 6;
n = 0;
for subj = subj_array
    n = n + 1;
    figure
    for i = 1:length(conds)
        rows = kinem.LateGroup.Subject==subj & strcmp(kinem.LateGroup.Type,conds{i});

        % Kinematics
        meanSway(n,i) = nanmean(kinem.LateGroup.StdSway(rows));
        meanCOMSway(n,i) = nanmean(kinem.LateGroup.StdCOMSway(rows));
        meanrFvertClav(n,i) = nanmean(kinem.LateGroup.rFvertClav(rows));
        meanrFlatCOM(n,i) = nanmean(kinem.LateGroup.rFlatCOM(rows));
        meanrFvertCOM(n,i) = nanmean(kinem.LateGroup.rFvertCOM(rows));
        meanrFlatClav(n,i) = nanmean(kinem.LateGroup.rFlatClav(rows));
        
        dist = kinem.LateGroup.Dist(rows);
    end
end

%% Force data
forces = load('HHI2017_LateStats_force_MW');
conds = {'Assist Ground','Assist Beam'};

%% Cycle through all subj's and plot parameter values for fitting force to int. pt. and clavicle. Do for Assist Beam only.
numrows = length(subj_array_force); numcols = 3;
plotind = 0;
xlab{1} = 'ML IP';
xlab{2} = 'Vert IP';
xlab{3} = 'ML Clav';

for subj = subj_array_force
    rows = forces.LateGroup.Subject==subj & strcmp(forces.LateGroup.Type,conds{2});
    % Int. Point - care ML and vert dir's
%     rxIP = forces.LateGroup.Rsqx(rows); 
%     mxIP = forces.LateGroup.mx(rows); 
%     bxIP = forces.LateGroup.bx(rows); 
%     kxIP = forces.LateGroup.kx(rows); 
%     rzIP = forces.LateGroup.Rsqz(rows); 
%     mzIP = forces.LateGroup.mz(rows); 
%     bzIP = forces.LateGroup.bz(rows); 
%     kzIP = forces.LateGroup.kz(rows); 
    % Clavicle - care ML dir only for sway
    rxClav = forces.LateGroup.Rsqx_clav(rows); 
    mxClav = forces.LateGroup.mx_clav(rows); 
    bxClav = forces.LateGroup.bx_clav(rows); 
    kxClav = forces.LateGroup.kx_clav(rows); 
    % Calculate percent of trials where param was sig.
%     mxIPsig = length(find(~isnan(mxIP)))/length(mxIP);
%     bxIPsig = length(find(~isnan(bxIP)))/length(bxIP);
%     kxIPsig = length(find(~isnan(kxIP)))/length(kxIP);
%     mzIPsig = length(find(~isnan(mzIP)))/length(mzIP);
%     bzIPsig = length(find(~isnan(bzIP)))/length(bzIP);
%     kzIPsig = length(find(~isnan(kzIP)))/length(kzIP);
    mxClavsig = length(find(~isnan(mxClav)))/length(mxClav);
    bxClavsig = length(find(~isnan(bxClav)))/length(bxClav);
    kxClavsig = length(find(~isnan(kxClav)))/length(kxClav);
    
%     % Plot all Rsq on same plot to compare which fits are best
%     plotind = plotind + 1;
%     subplot(numrows,numcols,plotind)
%     boxplot([rxIP rzIP rxClav],1:3),ylabel('R^2')
%     set(gca,'xticklabel',xlab);
%     titlename = sprintf('HHI%i',subj),title(titlename);
%     
%     plotind = plotind + 1;
%     subplot(numrows,numcols,plotind)
%     bar(1:3,[mxIPsig mzIPsig mxClavsig]),ylabel('% trials m sig'),hline(0,'k');
%     set(gca,'xticklabel',xlab);
%     
%     plotind = plotind + 1;
%     subplot(numrows,numcols,plotind)
%     bar(1:3,[bxIPsig bzIPsig bxClavsig]),ylabel('% trials b sig'),hline(0,'k');
%     set(gca,'xticklabel',xlab);
%     
%     plotind = plotind + 1;
%     subplot(numrows,numcols,plotind)
%     bar(1:3,[kxIPsig kzIPsig kxClavsig]),ylabel('% trials k sig'),hline(0,'k');
%     set(gca,'xticklabel',xlab);

    % Plot all mass values on same plot to get idea of range
    plotind = plotind + 1;
    subplot(numrows,numcols,plotind)
    boxplot([mxIP mzIP mxClav],1:3),ylabel('m (kg)'),hline(0,'k');
    set(gca,'xticklabel',xlab);
    titlename = sprintf('HHI%i',subj),title(titlename);
    
    % All damping values    
    plotind = plotind + 1;
    subplot(numrows,numcols,plotind)
    boxplot([bxIP bzIP bxClav],1:3),ylabel('b (N/(m/s))'),hline(0,'k');
    set(gca,'xticklabel',xlab);
    
    % All stiffness values    
    plotind = plotind + 1;
    subplot(numrows,numcols,plotind)
    boxplot([kxIP kzIP kxClav],1:3),ylabel('k (N/m)'),hline(0,'k');
    set(gca,'xticklabel',xlab);
    
end

%% Plot R^2 and each param value for each pair in x direction
subplot(4,1,1),boxplot(rxIP,subjArray),hold on,ylabel('R^2') 
subplot(4,1,2),boxplot(mxIP,subjArray),ylabel('m (kg)')
subplot(4,1,3),boxplot(bxIP,subjArray),ylabel('b (N/(m/s))')
subplot(4,1,4),boxplot(kxIP,subjArray),ylabel('k (N/m)'),xlabel('Partnership')

%% Plot % time where param value sig for each pair in x direction
subplot(3,1,1),bar(mxIPsig*100),ylim([0 100]),ylabel('% trials m sig'),box off, set(gca,'tickdir','out');
subplot(3,1,2),bar(bxIPsig*100),ylim([0 100]),ylabel('% trials b sig'),box off, set(gca,'tickdir','out');
subplot(3,1,3),bar(kxIPsig*100),ylim([0 100]),ylabel('% trials k sig'),xlabel('Partnership'),box off, set(gca,'tickdir','out');

%% Plot R^2 and each param value for each pair in z direction
subplot(4,1,1),boxplot(rzIP,subjArray),hold on,ylabel('R^2'),box off, set(gca,'tickdir','out');
subplot(4,1,2),boxplot(mzIP,subjArray),ylabel('m (kg)'),box off, set(gca,'tickdir','out'),hline(0,'k--')
subplot(4,1,3),boxplot(bzIP,subjArray),ylabel('b (N/(m/s))'),box off, set(gca,'tickdir','out'),hline(0,'k--')
subplot(4,1,4),boxplot(kzIP,subjArray),ylabel('k (N/m)'),xlabel,box off, set(gca,'tickdir','out'),('Partnership'),hline(0,'k--')

%% Plot % time where param value sig for each pair in z direction
subplot(3,1,1),bar(mzIPsig*100),ylim([0 100]),ylabel('% trials m sig'),box off, set(gca,'tickdir','out');
subplot(3,1,2),bar(bzIPsig*100),ylim([0 100]),ylabel('% trials b sig'),box off, set(gca,'tickdir','out');
subplot(3,1,3),bar(kzIPsig*100),ylim([0 100]),ylabel('% trials k sig'),xlabel('Partnership'),box off, set(gca,'tickdir','out');

%% Cycle through all subj's and plot mean and SD for xcorr and lag values
plotind = 0;
numrows = length(subj_array_force); numcols = 4;
for subj = subj_array_force
    clear xcorrX lagX xcorrZ lagZ
    for i = 1:length(conds)
        rows = forces.LateGroup.Subject==subj & strcmp(forces.LateGroup.Type,conds{1});
        % ML dir
        xcorrX(:,i) = forces.LateGroup.xcorrX(rows); 
        lagX(:,i) = forces.LateGroup.lagX(rows); 
        % Vert dir
        xcorrZ(:,i) = forces.LateGroup.xcorrZ(rows); 
        lagZ(:,i) = forces.LateGroup.lagZ(rows); 
    end
    
    plotind = plotind + 1;
    subplot(numrows,numcols,plotind)
    boxplot(xcorrX); set(gca,'xticklabel',conds); ylabel('ML F and disp xcorr norm'),...
        box off; set(gca,'tickdir','out');
        xlim([0.5 2.5]); ylim([-1 1]);
        titlename = sprintf('HHI%i',subj); title(titlename);
        
    plotind = plotind + 1;
    subplot(numrows,numcols,plotind)
    boxplot(lagX); set(gca,'xticklabel',conds); ylabel('ML F and disp lag (s)'),...
        box off; set(gca,'tickdir','out');
        xlim([0.5 2.5]); ylim([-1 1]);
        
    plotind = plotind + 1;
    subplot(numrows,numcols,plotind)
    boxplot(xcorrZ); set(gca,'xticklabel',conds); ylabel('Vert F and disp xcorr norm'),...
        box off; set(gca,'tickdir','out');
        xlim([0.5 2.5]); ylim([-1 1]);
        
    plotind = plotind + 1;
    subplot(numrows,numcols,plotind)
    boxplot(lagZ); set(gca,'xticklabel',conds); ylabel('Vert F and disp lag (s)'),...
        box off; set(gca,'tickdir','out');
        xlim([0.5 2.5]); ylim([-1 1]);
end

