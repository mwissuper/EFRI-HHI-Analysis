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
%       3. Calculate metrics measured for each
%       trial, and divide the statistics into tables
%       4. Save the processed data and the statistics table

%   Luke Drnach
%   December 7, 2018

%   Edited by M Wu 8/30/19-5/28/20 to look at additional variables/metrics

clear; close all; clc;

analyze_force = 1;
plotCheck = 0; % 1 = plot clav AP and ML pos to check distance completed on beam per trial. 2 = plot algo calc dist vs Luke's method. 3 = plot sway vs. beam midline

temp = cd;
sourcefolder = [temp,'\Data files'];
if analyze_force == 1
    subj_array = [3:5 8:13]; % Use 3-14 for kinem data, skip 6,7,14 for force data since expt error
else % More participants for kinem data only
    subj_array = 3:14; % 
end

%% Run analysis for work, power, sway per participant

for subj = 13%subj_array
    file = sprintf('HHI2017_%i.mat',subj);

    % Reorganize the data and calculate work and power - for an individual trial
    load([sourcefolder,'\',file]);
    TrialData = preprocessHHI2017(TrialData);                   % Reorganizes individual trial data
    % Manually correct marker labeling errors since don't know how to fix
    % in Nexus and then export to mat
    if subj == 3
        for n = [35 46 49]
            temp = TrialData(n).Markers.POB.LHEE;
            TrialData(n).Markers.POB.LHEE = TrialData(n).Markers.POB.LTOE;
            TrialData(n).Markers.POB.LTOE = temp;
        end
        n = 42;
        temp = TrialData(n).Markers.LHEE;
        TrialData(n).Markers.LHEE = TrialData(n).Markers.LTOE;
        TrialData(n).Markers.LTOE = temp;
    elseif subj == 4
        temp = TrialData(45).Markers.FH.frontleft;
        TrialData(45).Markers.FH.frontleft = TrialData(45).Markers.FH.backmiddle;
        TrialData(45).Markers.FH.backmiddle = TrialData(45).Markers.FH.frontright;
        TrialData(45).Markers.FH.frontright = TrialData(45).Markers.FH.frontmiddle;
        TrialData(45).Markers.FH.frontmiddle = temp;
    elseif subj == 14
        n = 31;
        temp = TrialData(n).Markers.POB.LHEE;
        TrialData(n).Markers.POB.LHEE = TrialData(n).Markers.POB.LTOE;
        TrialData(n).Markers.POB.LTOE = temp;
    end
    TrialData = mainWorkPowerAnalysisMW(TrialData,subj,plotCheck);               % Calculates work and power transfer
    if plotCheck == 1
        figname = sprintf('HHI%i check time window.fig',subj);
    elseif plotCheck == 2
        figname = sprintf('HHI%i check dist beam.fig',subj);
    end
%     saveas(gcf,figname,'fig');
%     close(gcf);
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
plotoption = 0; % Plot each component force with peaks
GroupStats = [];
for subj = subj_array 
    file = sprintf('HHI2017_%i.mat',subj) % file from current folder, not source data
    load(file); % Now use new file with processed data
    GroupStats = HHI2017GroupStats_MW(TrialData,GroupStats,subj);
end

%% Save the final results to current folder
% Save the Trial Data with Calculated Work and Power, and save the grouped
% statistics
fprintf('Saving Results\n');

if analyze_force == 1
    save('HHI2017_Stats_force_MW.mat','GroupStats');
else
    save('HHI2017_Stats_MW.mat','GroupStats');
end

% writetable(EarlyGroup,'HHI_Block1.xlsx');
% writetable(LateGroup,'HHI_Stats.xlsx');