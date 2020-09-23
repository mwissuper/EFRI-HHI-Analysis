function plotForceDistributions(TrialData)
%% PLOTFORCEDISTRIBUTIONS: Histogram plots of forces from a HHI Experiment
%
%   plotForceDistributions(TRIALDATA) makes histograms of the three
%   interaction forces recorded in TrialData.Results.Forces.
%   PLOTFORCEDISTRIBUTIONS collects the forces across all entries in
%   TrialData (i.e. across all listed trials) and makes histograms for the
%   Assist Beam and Assist Ground trial types.

%   Luke Drnach
%   December 6, 2018


%% Collect all of the forces from the TrialData structure
numTrials = length(TrialData);
allForces = cell(1,numTrials);
for n = 1:numTrials
    if ~isempty(TrialData(n).Results)
        allForces{n} = TrialData(n).Results.Forces;
    end
end
% Index which trials are Assist Beam and which are Assist Ground.
Info = [TrialData.Info];
Type = {Info.Condition};
AB_Trials = strcmpi(Type,'Assist Beam');
AG_Trials = strcmpi(Type,'Assist Ground');
%% Plot the force histograms
forceAB = cat(2,allForces{AB_Trials});
forceAG = cat(2,allForces{AG_Trials});
labels = {'Fx (N)','Fy (N)','Fz (N)'};
for n = 1:3
    figure()
    histogram(forceAB(n,:), 'FaceAlpha',0.5,'Normalization','pdf');
    hold on;
    histogram(forceAG(n,:), 'FaceAlpha',0.5,'Normalization','pdf');
    xlabel(labels{n});
    set(get(gcf,'Children'), 'FontName','Garamond','FontSize',12,'FontWeight','bold');
    legend('Beam','Ground','location','northwest');
    legend boxoff;
end
end

