function plotHHIStatsAggregated(stats)
%% PLOTHHISTATSAGGREGATED: Plots Statistics from HHI Experiment aggregated by type
%
%   plotHHIStatsAggregated(STATS) takes a table of statistics STATS from
%   the Human-Human Interaction experiment and produces boxplots of Peak
%   work, power, and force values in assisted beam and ground walking.
%   PLOTHHISTATSAGGREGATED also makes boxplots for Distance Traveled, Average
%   Speed, and Sway for all four trial types.
%
%   PLOTHHISTATSAGGREGATED aggregates the values across all listed
%   participants. To get plots for a single participant, the input table
%   STATS must contain only one participant pair.

%   Luke Drnach
%   December 6, 2018

% First get the indices for each trial type
types = stats.Type(:);
idxAB = strcmpi(types,'Assist Beam');
idxAG = strcmpi(types,'Assist Ground');
idxSB = strcmpi(types,'Solo Beam');
idxSG = strcmpi(types,'Solo Ground');

%% Aggregate Data plots
% Preliminary Calculations: Number of Assist Beam and Assist Ground
% Trials
numAB = sum(idxAB);
numAG = sum(idxAG);
maxNum = max([numAB, numAG]);
%% Work and Power Boxplots
negfields = {'PeakNegPower','PeakNegWork','PeakZCompression','PeakNegFx','PeakNegFy'};
posfields = {'PeakPosPower','PeakPosWork','PeakZTension','PeakPosFx','PeakPosFy'};
labels = {'Peak Power (W)', 'Peak Work (J)','Peak Vertical Force (N)','Peak Lateral Force (N)','Peak Forward Force (N)'};
for n = 1:length(labels)
    figure();
    % Positive Values
    subplot(1,2,1)
    data = NaN(maxNum,2);
    data(1:numAG,1) = stats.(posfields{n})(idxAG);
    data(1:numAB,2) = stats.(posfields{n})(idxAB);
    boxplot(data);
    ylabel(['Positive ',labels{n}]);
    %Negative Values
    subplot(1,2,2)
    data = NaN(maxNum,2);
    data(1:numAG,1) = stats.(negfields{n})(idxAG);
    data(1:numAB,2) = stats.(negfields{n})(idxAB);
    boxplot(data);
    ylabel(['Negative ',labels{n}]);
    set(get(gcf,'Children'), 'FontName','Garamond','FontSize',12,'FontWeight','bold',...
        'XTickLabel',{'Ground','Beam'},'XTickLabelRotation',45)
end

%% BeamWalking Performance Boxplots
% Plot Distance Traveled, Average Speed, and Sway as boxplots for all four
% trial types (not including Assist Solo)
fields = {'Dist','AvgSpeed','StdSway'};
labels = {'Distance Traveled (m)','Average Speed (m/s)','Sway (mm)'};
for n = 1:length(fields)
   figure();
   SG = stats.(fields{n})(idxSG);
   AG = stats.(fields{n})(idxAG);
   SB = stats.(fields{n})(idxSB);
   AB = stats.(fields{n})(idxAB);
   boxplot([SG,AG,SB,AB]);
   ylabel(labels{n});
   set(gca,'XTickLabel',{'SG','AG','SB','AB'},'FontName','Garamond','FontSize',12,'FontWeight','bold');
end
end