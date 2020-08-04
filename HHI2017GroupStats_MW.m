function [GroupStats] = HHI2017GroupStats_MW(TrialData,GroupStats,subj)
%% Calculates the statistics using all trials of the experiment
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

%% Calculate the statistics
% Get the statistics for the first block, and append them to the
% existing statistics if already did first subj since now there are more
% columns of data won't concatenate, also don't need previous data repeated
plotoption = 0;
Stat = getPointStatsHHI_MW(TrialData(idx),plotoption); % don't plot
if subj > 3
    GroupStats = [GroupStats; Stat];
else % Special case of first subject, need to set up array size
    GroupStats = Stat;
end

