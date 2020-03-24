function plotHHIStatsBySubject(stats)
    %% PLOTHHISTATSBYSUBJECT: Makes plots of HHI statistics grouped by subject
    %
    %   plotHHIStatsBySubject(STATS) plots peak power, work, and forces, as
    %   well as performance measures, in different figures.
    %   plotHHIStatsBySubject plots the mean and standard deviation of each
    %   measure, color-coded by subject number
    
    %   Luke Drnach
    %   January 8, 2018
    %
    %   Updated and refactored on December 6, 2018. Updated to use tables
    %   as the main format for handling point statistics. Refactored for
    %   code clarity.
    
    %% List the plots to be made HERE (EDIT THIS PART FOR NEW PLOTS)
    % For comparisons involving only the ASSISTED trials (assist Ground
    % and assist beam), include the variable in ASSISTFIELDS and
    % provide and appropriate label in ASSISTLABELS
    assistFields = {'PeakNegPower','PeakNegWork','PeakZCompression','PeakPosFx',...
        'PeakPosFy'};
    assistLabels = {'Max Negative Power (W)','Max Negative Work (J)',...
        'Max Compressive Force (N)', 'Max Lateral Force (N)',...
        'Max Forward Force (N)'};
    
    % For comparisons involving ALL trial conditions (except assist
    % solo), include the variable for comparison in ALLFIELDS, and
    % provide an appropriate label in ALLLABELS
    allFields = {'Dist','AvgSpeed','StdSway'};
    allLabels = {'Distance Traveled (m)','Average Speed (m/s)','Sway (mm)'};
    
    %% Set-up for plotting
    % Get the number of unique subjects
    subjects = unique(stats.Subject(:));
    lgdstr = strsplit(num2str(subjects'));   % Subject numbers for the legend
    numSubjects = length(lgdstr);                           % Number of subjects
    
    % Set up the plotting
    nAssist = length(assistFields);
    numPlots = nAssist + length(allFields);
    F(1:numPlots+1) = figure();
    for f = 2:numPlots+1
        F(f) = figure();
    end
    MSize = 7; % Marker Size
    
    % Concatenate the fields and labels into one cell array
    fields = [assistFields,allFields];
    labels = [assistLabels,allLabels];
    
    % Get the indices for each trial type
    types = stats.Type(:);
    idxAB = strcmpi(types,'Assist Beam');
    idxAG = strcmpi(types,'Assist Ground');
    idxSB = strcmpi(types,'Solo Beam');
    idxSG = strcmpi(types,'Solo Ground');
    %% Plot the assisted (2-point comparison) trials
    % These plots are just for assisted trials (2-point comparison)
    for k = 1:nAssist
        figure(F(k));
        for m = 1:numSubjects
            mAG = mean(stats.(fields{k})(and(idxAG,stats.Subject==subjects(m))),'omitnan');
            mAB = mean(stats.(fields{k})(and(idxAB,stats.Subject==subjects(m))),'omitnan');
            sAG = std(stats.(fields{k})(and(idxAG,stats.Subject==subjects(m))));
            sAB = std(stats.(fields{k})(and(idxAB,stats.Subject==subjects(m))));
            errorbar([1,2],[mAG,mAB],[sAG,sAB],'-o','MarkerSize',MSize,'MarkerFaceColor','auto');
            hold on;
        end
        set(gca,'XTick',1:2,'XTickLabel',{'AG','AB'},'FontName','Garamond','FontSize',14,'FontWeight','bold');
        ylabel(labels{k});
        xlim([0,3]);
    end
    %% Plot the 4-way comparisons
    for k = nAssist+1:numPlots
        figure(F(k));
        for m = 1:numSubjects
            % Solo Ground
            mSG = mean(stats.(fields{k})(and(idxSG,stats.Subject==subjects(m))));
            sSG = std(stats.(fields{k})(and(idxSG,stats.Subject==subjects(m))));
            % Assist Ground
            mAG = mean(stats.(fields{k})(and(idxAG,stats.Subject==subjects(m))));
            sAG = std(stats.(fields{k})(and(idxAG,stats.Subject==subjects(m))));
            % Solo Beam
            mSB = mean(stats.(fields{k})(and(idxSB,stats.Subject==subjects(m))));
            sSB = std(stats.(fields{k})(and(idxSB,stats.Subject==subjects(m))));
            % Assist Beam
            mAB = mean(stats.(fields{k})(and(idxAB,stats.Subject==subjects(m))));
            sAB = std(stats.(fields{k})(and(idxAB,stats.Subject==subjects(m))));
            % PLOT
            errorbar([1,2,3,4],[mSG,mAG,mSB,mAB],[sSG,sAG,sSB,sAB],'-o','MarkerSize',MSize,'MarkerFaceColor','auto');
            hold on;
        end
        set(gca,'XTick',1:4,'XTickLabel',{'SG','AG','SB','AB'},'FontName','Garamond','FontSize',14,'FontWeight','bold');
        ylabel(labels{k});
        xlim([0,5])
    end
    %% Make a figure just for a subject color code
    figure(F(end))
    for n = 1:numSubjects
        plot(0,0,'.');
        hold on;
    end
    legend(lgdstr{:})
    legend boxoff;
end
