%% PIPELINEHHI: Analysis pipeline for Human-Human Interaction Experiment (2017)
%
%   PIPELINEHHI is a script that runs a sequence of functions for analyzing
%   the data from the 2017 HHI experiment. PIPELINEHHI loads the data
%   (processed by concatData) and the previous statistics. The steps of the
%   pipeline are:
%       1. Preprocess to remove unknown marker fields and subject-specific
%       identifiers in the Marker Data
%       2. Calculate rectified forces, work, power, and performance
%       measured for each trial
%       3. Replace the repeated trials in the dataset
%       4. Calculate the peak values and other point measured for each
%       trial, and divide the statistics into tables for the first block
%       and the remaining trials.
%       5. Save the processed data and the statistics table
%       6. Write the statistics to an EXCEL file.

%   Luke Drnach
%   December 7, 2018

%   Edited by M Wu 8/30/19-1/31/20 to look at additional variables/metrics

clear; close all; clc;

analyze_force = 1;

temp = cd;
sourcefolder = [temp,'\Data files'];
if analyze_force == 1
    subj_array = [3:5 8:13]; % Use 3-14 for kinem data, skip 6,7,14 for force data since expt error
else % More participants for kinem data only
    subj_array = 3:14; % 
end

%% Run analysis for work, power, sway per participant

for subj = 3%subj_array
    file = sprintf('HHI2017_%i.mat',subj);

    % Reorganize the data and calculate work and power - for an individual trial
    load([sourcefolder,'\',file]);
    TrialData = preprocessHHI2017(TrialData);                   % Reorganizes individual trial data
    TrialData = mainWorkPowerAnalysisMW(TrialData,subj);               % Calculates work and power transfer
    % Some don't have NOTES. Check for NOTES
    for n = 1:length(TrialData)
        if ~isfield(TrialData(n).Info,'Notes')
            TrialData(n).Info.Notes = ' ';
        end
    end
    filtTrialData = filterRepeatTrials(TrialData);
    save(file,'TrialData'); % Same file names as raw trials in sourcefolder but save to current folder
end

%% Loop through all metrics for individual participants and concatenate for stats and plots
% Do only for subj's with good force data
EarlyGroup = []; LateGroup = [];
% plotoption = 1; % Plot vector power with peaks
% plotoption = 2; % Plot each component power with peaks
plotoption = 0; % Plot each component force with peaks
for subj = subj_array 
    file = sprintf('HHI2017_%i.mat',subj) % file from current folder, not source data
    load(file); % Now use new file with processed data
    [EarlyGroup, LateGroup] = HHI2017EarlyLateGroupStats_MW(TrialData,EarlyGroup,LateGroup,subj,plotoption);
end

%% Save the final results to current folder
% Save the Trial Data with Calculated Work and Power, and save the grouped
% statistics
fprintf('Saving Results\n');

if analyze_force == 1
    save('HHI2017_EarlyStats_force_MW.mat','EarlyGroup');
    save('HHI2017_LateStats_force_MW.mat','LateGroup');
else
    save('HHI2017_EarlyStats_MW.mat','EarlyGroup');
    save('HHI2017_LateStats_MW.mat','LateGroup');
end

% %% Export the statistics to EXCEL
% writetable(EarlyGroup,'HHI_Block1.xlsx');
% writetable(LateGroup,'HHI_Stats.xlsx');