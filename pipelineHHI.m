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

%   Edited by M Wu 8/30/19 to look at additional variables/metrics

%% Current file to run
file = 'HHI2017_13.mat';

%% Paths for group statistics and saving
sourcefolder = 'C:\Users\mwu40\Documents\EFRI-HHI\EFRI-HHI Data';
LateFile = 'HHI2017_LateStats.mat';
EarlyFile = 'HHI2017_EarlyStats.mat';

%% Reorganize the data and calculate work and power - for an individual trial
load([sourcefolder,'\',file]);
TrialData = preprocessHHI2017(TrialData);                   % Reorganizes individual trial data
TrialData = mainWorkPowerAnalysisMW(TrialData);               % Calculates work and power transfer
% Some don't have NOTES. Check for NOTES
for n = 1:length(TrialData)
    if ~isfield(TrialData(n).Info,'Notes')
        TrialData(n).Info.Notes = ' ';
    end
end
filtTrialData = filterRepeatTrials(TrialData);
%% Calculate the Group Statistics
fprintf('Processing Group Statistics\n');
load([sourcefolder,'\',LateFile]);
load([sourcefolder,'\',EarlyFile]);
[EarlyGroup, LateGroup] = HHI2017EarlyLateGroupStats(filtTrialData,EarlyGroup,LateGroup);

%% Save the final results
% Save the Trial Data with Calculated Work and Power, and save the grouped
% statistics
fprintf('Saving Results\n');
save([sourcefolder,'\',file],'TrialData');
save([sourcefolder,'\',EarlyFile],'EarlyGroup');
save([sourcefolder,'\',LateFile],'LateGroup');

% %% Export the statistics to EXCEL
% writetable(EarlyGroup,'HHI_Block1.xlsx');
% writetable(LateGroup,'HHI_Stats.xlsx');