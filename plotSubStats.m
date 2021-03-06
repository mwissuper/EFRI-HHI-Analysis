% plot individual subj's stats data using concatenated tables from
% LateGroup
clear all; clc; close all;

subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)
plotPosNeg = 1;
if plotPosNeg == 1
    xlims = [-12 12];
    numcols = 3; % Look all dir's force
    numrows = 3;
else
    xlims = [0 15];
    numcols = 3; % Only look Fvert
    numrows = 3;
end

%% Kinematics
kinem = load('HHI2017_LateStats_MW');
conds = {'Solo Ground','Assist Ground','Solo Beam','Assist Beam'};

n = 0;
for subj = subj_array
    n = n + 1;
    for i = 1:length(conds)
        rows = kinem.LateGroup.Subject==subj & strcmp(kinem.LateGroup.Type,conds{i});

        % Kinematics
        StdSway(n,i) = nanmean(kinem.LateGroup.StdSway(rows));
        StdCOMSway(n,i) = nanmean(kinem.LateGroup.StdCOMSway(rows));
        rFvertClav(n,i) = nanmean(kinem.LateGroup.rFvertClav(rows));
        rFlatCOM(n,i) = nanmean(kinem.LateGroup.rFlatCOM(rows));
        rFvertCOM(n,i) = nanmean(kinem.LateGroup.rFvertCOM(rows));
        rFlatClav(n,i) = nanmean(kinem.LateGroup.rFlatClav(rows));

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
    rxIP = forces.LateGroup.Rsqx(rows); 
    mxIP = forces.LateGroup.mx(rows); 
    bxIP = forces.LateGroup.bx(rows); 
    kxIP = forces.LateGroup.kx(rows); 
    rzIP = forces.LateGroup.Rsqz(rows); 
    mzIP = forces.LateGroup.mz(rows); 
    bzIP = forces.LateGroup.bz(rows); 
    kzIP = forces.LateGroup.kz(rows); 
    % Clavicle - care ML dir only for sway
    rxClav = forces.LateGroup.Rsqx_clav(rows); 
    mxClav = forces.LateGroup.mx_clav(rows); 
    bxClav = forces.LateGroup.bx_clav(rows); 
    kxClav = forces.LateGroup.kx_clav(rows); 
    % Calculate percent of trials where param was sig.
    mxIPsig = length(find(~isnan(mxIP)))/length(mxIP);
    bxIPsig = length(find(~isnan(bxIP)))/length(bxIP);
    kxIPsig = length(find(~isnan(kxIP)))/length(kxIP);
    mzIPsig = length(find(~isnan(mzIP)))/length(mzIP);
    bzIPsig = length(find(~isnan(bzIP)))/length(bzIP);
    kzIPsig = length(find(~isnan(kzIP)))/length(kzIP);
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

%% Cycle through all subj's and concatenate parameter values for fitting force to int. pt. Do for Assist Beam only. param value vs. subj
numrows = length(subj_array_force);
subj_ind = 0;
subjArray = [];
rxIP = []; mxIP = []; bxIP = []; kxIP = [];
mxIPsig = []; bxIPsig = []; kxIPsig = [];
rzIP = []; mzIP = []; bzIP = []; kzIP = [];
mzIPsig = []; bzIPsig = []; kzIPsig = [];
% Concatenate per subj
for subj = subj_array_force
    subj_ind = subj_ind + 1;
    rows = forces.LateGroup.Subject==subj & strcmp(forces.LateGroup.Type,conds{2});
    a = find(rows~=0);
    subjArray(end+1:end+length(a),1) = subj_ind;
    % Int. Point - care ML and vert dir's
    rxIP = [rxIP; forces.LateGroup.Rsqx(rows)]; 
    mxIP = [mxIP; forces.LateGroup.mx(rows)]; 
    bxIP = [bxIP; forces.LateGroup.bx(rows)]; 
    kxIP = [kxIP; forces.LateGroup.kx(rows)]; 
    
    rzIP = [rzIP; forces.LateGroup.Rsqz(rows)]; 
    mzIP = [mzIP; forces.LateGroup.mz(rows)]; 
    bzIP = [bzIP; forces.LateGroup.bz(rows)]; 
    kzIP = [kzIP; forces.LateGroup.kz(rows)];
    
    % Percent rials where param sig
    % Calculate percent of trials where param was sig.
    mxIPsig = [mxIPsig; length(find(~isnan(forces.LateGroup.mx(rows))))/length(forces.LateGroup.mx(rows))];
    bxIPsig = [bxIPsig; length(find(~isnan(forces.LateGroup.bx(rows))))/length(forces.LateGroup.bx(rows))];
    kxIPsig = [kxIPsig; length(find(~isnan(forces.LateGroup.kx(rows))))/length(forces.LateGroup.kx(rows))];
    mzIPsig = [mzIPsig; length(find(~isnan(forces.LateGroup.mz(rows))))/length(forces.LateGroup.mz(rows))];
    bzIPsig = [bzIPsig; length(find(~isnan(forces.LateGroup.bz(rows))))/length(forces.LateGroup.bz(rows))];
    kzIPsig = [kzIPsig; length(find(~isnan(forces.LateGroup.kz(rows))))/length(forces.LateGroup.kz(rows))];
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

%% Cycle through all subj's and plot histogram for Assist Beam Fx or Fz (just comment part of code)
% Sturge's rule says choose # of bins by 1 + 3.322*log(n), where n is
% number of observations. For each subject, 8 trials, so about 8 bins. Bin
% width of 0.5 accomplishes this.
plotind = 0;
subjind = 0;
for subj = subj_array_force
    clear Fx Fy Fz
    subjind = subjind + 1;
    rows = forces.LateGroup.Subject==subj & strcmp(forces.LateGroup.Type,conds{2});
%     length(find(rows~=0))
    if plotPosNeg == 1
        % Combine pos and neg means into one data set for histogram
        Fx = [forces.LateGroup.meanPosFx(rows); forces.LateGroup.meanNegFx(rows)]; 
        Fy = [forces.LateGroup.meanPosFy(rows); forces.LateGroup.meanNegFy(rows)]; 
        Fz = [forces.LateGroup.meanPosFz(rows); forces.LateGroup.meanNegFz(rows)]; 
    else
        % Force mag (unsigned)
        Fx = forces.LateGroup.meanFx(rows); 
        Fy = forces.LateGroup.meanFy(rows);
        Fz = forces.LateGroup.meanFz(rows);
    end

    %% Force histogram. Matlab's auto param's chosen to reveal shape of distribution. If comparing subj's, may need to fix bin width...

    % Fx
    plotind = plotind + 1;
    subplot(numrows,numcols,plotind)
    histogram(Fx,'binwidth',0.5,'Normalization','probability'); box off; set(gca,'tickdir','out');
    xlim(xlims); ylabel('Probability (norm cts)');
    ylim([0 0.3]);
%     s = [859.0000  422.3333  574.6667  716.6667];
%     set(gcf,'outerposition',s);
    if plotPosNeg == 1
        xlabel('Mean mediolateral force (N)');
    else
        xlabel('Mean Fmag ML (N)');
    end
    titlename = sprintf('Partnership %i',subjind); title(titlename);

%     % Fz
%     plotind = plotind + 1;
%     subplot(numrows,numcols,plotind)
%     histogram(Fz,'binwidth',1,'Normalization','probability'); box off; set(gca,'tickdir','out');
%     xlim(xlims); ylabel('Probability (norm cts)');
%     ylim([0 0.8]);
%     s = [859.0000  422.3333  574.6667  716.6667];
%     set(gcf,'outerposition',s);
%     xlabel('Mean Fmag Vert (N)');
    
end

%% Cycle through all subj's and plot mean and SD for xcorr and lag values
forces = load('HHI2017_LateStats_force_MW');
conds = {'Assist Ground','Assist Beam'};
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

