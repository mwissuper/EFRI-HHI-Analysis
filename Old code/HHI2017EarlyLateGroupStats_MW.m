function [EarlyGroup,LateGroup] = HHI2017EarlyLateGroupStats_MW(TrialData,EarlyGroup,LateGroup,subj,plotoption)
%% HHI2017EarlyLateStats: Calculates the statistics for the early and late phases of the experiment
%
%   Early: The first two trials of each type
%   Late: The remaining trials of each type.


%   Luke Drnach
%   May 14, 2018

%   Updated December 7, 2018 to work with refactored HHI code

% Edited MW 9/5/19 to include new metrics and not to append unless we've
% done one participant

Info = [TrialData.Info];
Conditions = {Info.Condition};

idx = 1:length(Conditions);

% Find the indices of the first two of each trial type (will almost always
% be 1-10, but this will account for errors in running the experiment).
EarlyIdx = zeros(1,10);
C = {'Assist Beam','Assist Ground','Assist Solo','Solo Beam','Solo Ground'};
k = 1;
for n = 1:5
    ii = find(strcmp(Conditions,C{n}));
    EarlyIdx(k:k+1) = ii(1:2);
    k = k+2;
end
EarlyIdx = sort(EarlyIdx);
LateIdx = idx;
LateIdx(EarlyIdx) = [];

%% Calculate the statistics
% Get the statistics for the first block, and append them to the
% existing statistics if already did first subj since now there are more
% columns of data won't concatenate, also don't need previous data repeated
% EarlyStat = getPointStatsHHI_MW(TrialData(EarlyIdx),0); % Don't plot the early trials
% Get the statistics for the remaining blocks, and append them to the
% existing statistics if we did first subj already.
LateStat = getPointStatsHHI_MW(TrialData(LateIdx),plotoption);
s = [-7          33        2576        1416];
if plotoption == 1
    figname = sprintf('HHI%i_PowerVec_LateTrials.fig',subj);
    set(gcf,'outerposition',s);
    saveas(gcf,figname,'fig');
    close gcf;
elseif plotoption == 2
    figname = sprintf('HHI%i_PowerComp_LateTrials.fig',subj);
    set(gcf,'outerposition',s);
    saveas(gcf,figname,'fig');
    close gcf;
elseif plotoption == 3
    figname = sprintf('HHI%i_ForceComp_LateTrials.fig',subj);
    set(gcf,'outerposition',s);
    saveas(gcf,figname,'fig');
    close gcf;
end
if subj > 3
%     EarlyGroup = [EarlyGroup; EarlyStat];
    LateGroup = [LateGroup; LateStat];
else % Special case of first subject, need to set up array size
%     EarlyGroup = EarlyStat;
    LateGroup = LateStat;
end

end
