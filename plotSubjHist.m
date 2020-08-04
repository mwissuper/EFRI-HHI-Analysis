% Plot force and power histograms for an individual subject. Also
% scatterplot mean force/power metrics for subj vs. solo beam distance

clear; clc; close all;
subj_array = 3:14;
subj_array_force = [3:5 8:13];

colors = hsv(12);
colors(3,:) = [0.85,0.33,0.10]; % replace yellow
colors(5,:) = [0.47,0.67,0.19]; % replace a green
colors(6,:) = [0.5,0.5,0.5]; 
colors(8,:) = [0.00,0.50,1.00];
colors(12,:) = [0.49,0.18,0.56];

plotHist = 2; % 1 = force, 2 = power vs. another metric for scatterplots
plotMetric = 2; % 1 = solo beam distance, 2 = change in sway var
numcols = 3; numrows = 3;

plotind = 0;

%% Order subj's by performance in solo beam-walking

% Concatenate performance metric data
kinem = load('HHI2017_Stats_MW.mat');
% Find means for each trial type for each subj
n = 0;
subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)
conds_perf = {'Solo Beam','Assist Beam'};

n = 0;
for subj = subj_array % HHI12 had close to mean sway, good example subj
    n = n + 1;
    for i = 1:length(conds_perf)
        rows = kinem.GroupStats.Subject==subj & strcmp(kinem.GroupStats.Type,conds_perf{i});
        Dist(n,i) = nanmean(kinem.GroupStats.Dist(rows)); % Total distance traveled on beam
        StdSway(n,i) = nanmean(kinem.GroupStats.StdSway(rows)); % Sway variability
    end
end

[c,ia,ib] = intersect(subj_array,subj_array_force);
soloDistF = Dist(ia,1); % in ascending order for subj's with F data
dSwayV = StdSway(:,2) - StdSway(:,1);
dSwayVF = dSwayV(ia);
if plotMetric == 1 % sort histograms by solo beam dist
    [soloDistFSort,iSort] = sort(soloDistF); % sort worst to best
else % sort by change in sway var
    [swayVFSort,iSort] = sort(dSwayVF); % sort worst to best
end
kinem = [];

%% Concatenate all force data for each subj
n = 0;
for subj = subj_array_force
    n = n + 1;
    filename = sprintf('HHI2017_%i.mat',subj);
    load(filename);
    Fxtot = []; Pxtot = [];
    for i = 1:length(TrialData)
        if strcmp(TrialData(i).Info.Condition,'Assist Beam')
            Fxtot = [Fxtot; TrialData(i).Results.Forces(:,1)]; % Concatenate data together into one big array  
            Pxtot = [Pxtot; TrialData(i).Results.IntPower(:,1)]; 
        end
    end
    Fx{n} = Fxtot;
    Px{n} = Pxtot;
    % Keep track of which subj has smallest range
    rangeFx(n) = max(Fx{n}) - min(Fx{n});
    rangePx(n) = max(Px{n}) - min(Px{n});
    len(n) = length(Fx{n});
end

%% Choose hist parameters. Matlab's auto param's chosen to reveal shape of 
% distribution, but looks too coarse. Sturges rule specifies k bins
% from k = log2(n) + 1 where n = length. Keep track of bin width calculated
% for each subj based on Sturges rule.
n = 0;
for subj = subj_array_force
    n = n + 1;
    k(n) = ceil(log2(len(n))+1); % num bins
    if plotHist == 1
        h = histogram(Fx{n},k(n),'binmethod','Sturges'); % Make histogram just to get binwidth
    else
        h = histogram(Px{n},k(n),'binmethod','Sturges'); % Make histogram just to get binwidth
    end
    bw(n) = get(h,'binwidth');
    close all;
end

%% Plot histograms using median binwidth calc from Sturges rule
% Output metrics rangeHist, pkF, and P0 are sorted in order corresponding
% to subj's sorted by solo distance on beam
n = 0;
for subj = subj_array_force(iSort)
    n = n + 1;
    subplot(numrows,numcols,n)
    % Find color
    ind = find(subj_array == subj,1,'first'); % Use same color as kinem data
    if plotHist == 1
        h = histogram(Fx{n},'binwidth',median(bw),'Normalization','probability','facecolor',colors(ind,:),'edgecolor',colors(ind,:)); 
        rangeF(n) = h.BinLimits(2) - h.BinLimits(1);
        [maxP,i] = max(h.Values);
        pkF(n) = mean([h.BinEdges(i) h.BinEdges(i+1)]);
    else
        h = histogram(Px{n},'binwidth',median(bw),'Normalization','probability','facecolor',colors(ind,:),'edgecolor',colors(ind,:)); 
        ind = find(h.BinEdges == 0,1,'first');
        P0(n) = h.Values(ind);
        rangeP(n) = h.BinLimits(2) - h.BinLimits(1);
    end
    box off; set(gca,'tickdir','out'); % When use 'binmethod' of 'sturges' I get different num bins for each patnership! More consistent to specify bin width k
    titlename = sprintf('HHI %i',subj); title(titlename);
%     xlim([-39 39]),ylim([0 0.5]); 
    if plotHist == 1
        xlim([-15 15]);ylim([0 .4]);
    else
        xlim([-3 3]);ylim([0 .5]);
    end
    vline(0,'k--');
    if n == 8
        if plotHist == 1
            xlabel('Mediolateral force (N)');
        else
            xlabel('Mediolateral power (W)');
        end
    elseif n == 4
        ylabel('Probability');
    end
end

if plotMetric == 1
    x = soloDistF;
else
    x = dSwayVF;
end
        
if plotHist == 2 % Power metrics
    %% Plot probability of zero power vs. solo beam distance to see if people with better balance are more likely to not use energy exchange
    % histogram metrics are already in order of subj's solo distance
    % ability, so need to plot it this way. Need to make special legend
       
    figure
    subplot(1,2,1)
    for i = 1:length(iSort)
        subj = subj_array_force(iSort(i));
        ind = find(subj_array_force == subj,1,'first');
        indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
        plot(x(ind),rangeP(i),'.','markersize',14,'color',colors(indC,:)),hold on
    end
%     temp = {'3','4','5','8','9','10','11','12','13'};
%     legend(temp{iSort})
    if plotMetric == 1
        xlabel('Solo beam dist (m)');
    else
        xlabel('Change Sway Var (m)');
    end
    ylabel('P range (W)'),box off

    subplot(1,2,2)
    for i = 1:length(iSort)
        subj = subj_array_force(iSort(i));
        ind = find(subj_array_force == subj,1,'first');
        indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
        plot(x(ind),P0(i),'.','markersize',14,'color',colors(indC,:)),hold on
    end
    if plotMetric == 1
        xlabel('Solo beam dist (m)');
    else
        xlabel('Change Sway Var (m)');
    end
    ylabel('Prob. 0 power'),box off
else
    %% Force metrics
    figure
    subplot(1,2,1)
    for i = 1:length(subj_array_force)
        subj = subj_array_force(iSort(i));
        ind = find(subj_array_force == subj,1,'first');
        indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
        plot(x(ind),rangeF(i),'.','markersize',14,'color',colors(indC,:)),hold on;
    end
    if plotMetric == 1
        xlabel('Solo beam dist (m)');
    else
        xlabel('Change Sway Var (m)');
    end
    ylabel('F range (N)'),box off

    subplot(1,2,2)
    for i = 1:length(iSort)
        subj = subj_array_force(iSort(i));
        ind = find(subj_array_force == subj,1,'first');
        indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
        plot(x(ind),pkF(i),'.','markersize',14,'color',colors(indC,:)),hold on;
    end
    if plotMetric == 1
        xlabel('Solo beam dist (m)');
    else
        xlabel('Change Sway Var (m)');
    end
    ylabel('Most probable F (N)'),box off
end
set(gcf,'outerposition',[516   376   576   311]);