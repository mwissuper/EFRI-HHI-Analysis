% Plot force and power histograms for each partnership

clear;
clc; close all;
subj_array = 3:14;
subj_array_force = [3:5 8:13];

colors = hsv(12);
colors(3,:) = [0.85,0.33,0.10]; % replace yellow
colors(5,:) = [0.47,0.67,0.19]; % replace a green
colors(6,:) = [0.5,0.5,0.5]; 
colors(8,:) = [0.00,0.50,1.00];
colors(12,:) = [0.49,0.18,0.56];

var = 2; % 1 = force, 2 = power vs. another metric for scatterplots
numcols = 3; numrows = 3;

plotind = 0;

%% Order subj's by reduction in sway from solo to partnered beam-walking

% Concatenate performance metric data
kinem = load('HHI2017_Stats_MW.mat');
% Find means for each trial type for each subj
conds_perf = {'Solo Beam','Assist Beam'};

n = 0;
for subj = subj_array 
    n = n + 1;
    for i = 1:length(conds_perf)
        rows = kinem.GroupStats.Subject==subj & strcmp(kinem.GroupStats.Type,conds_perf{i});
        Dist(n,i) = nanmean(kinem.GroupStats.Dist(rows)); % Total distance traveled on beam
        StdSway(n,i) = nanmean(kinem.GroupStats.StdSway(rows)); % Sway variability
    end
end

[c,ia,ib] = intersect(subj_array,subj_array_force);

dSwayV = StdSway(:,2) - StdSway(:,1);
dSwayVF = dSwayV(ia);
[swayVFSort,iSort] = sort(dSwayVF); % sort worst to best

kinem = [];

%% Concatenate all force data for each subj
n = 0;
for subj = subj_array_force
    n = n + 1;
    filename = sprintf('HHI2017_%i.mat',subj);
    load(filename);
    Fxtot = []; Pxtot = []; tempFx = []; tempPx = [];
    for i = 1:length(TrialData)
        if strcmp(TrialData(i).Info.Condition,'Assist Beam')
            Fxtot = [Fxtot; TrialData(i).Results.Forces(:,1)]; % Concatenate data together into one big array  
            tempFx = [tempFx; nanmean(TrialData(i).Results.Forces(:,1))]; % Check that subj is correct by match the mean here to color and mean in mean F figure
            Pxtot = [Pxtot; TrialData(i).Results.IntPower(:,1)]; 
            tempPx = [tempPx; nanmean(TrialData(i).Results.IntPower(:,1))]; % Check that subj is correct by match the mean here to color and mean in mean P figure
        end
    end
    meanFx(n) = nanmean(tempFx);
    Fx{n} = Fxtot;
    meanPx(n) = nanmean(tempPx);
    Px{n} = Pxtot;
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
    if var == 1
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
    if var == 1
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
    if var == 1
        xlim([-15 15]);ylim([0 .4]);
    else
        xlim([-3 3]);ylim([0 .5]);
    end
    vline(0,'k--');
    if n == 8
        if var == 1
            xlabel('Mediolateral force (N)');
        else
            xlabel('Mediolateral power (W)');
        end
    elseif n == 4
        ylabel('Probability');
    end
end
