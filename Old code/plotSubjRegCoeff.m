% Pull out particular subject and look at regression fits to force in
% detail by plotting CI's on parameter values for all Assist Beam trials

clear all; clc; close all;

subj = 12;
file = sprintf('HHI2017_%i.mat',subj);
load(file);

% Pull out Assist Beam trials z dir only
coeff{1} = 'm';
coeff{2} = 'b';
coeff{3} = 'k';

numrows = 2; numcols = 5;
plotind = 0;
for i = 1:length(TrialData)
    if strcmp(TrialData(i).Info.Condition,'Assist Beam') 
        plotind = plotind + 1;
        subplot(numrows,numcols,plotind)
        errorbar(1:3,TrialData(i).Results.cintx(1:3,1),TrialData(i).Results.cintx(1:3,2))
        xlim([0.5 3.5]),set(gca,'xtick',1:3,'xticklabel',coeff);
        hline(0,'k--');
        if subj == 8
            ylim([-100 25]);
        else
            ylim([-80 5]);
        end
        
        if plotind == 1
            titlename = sprintf('%s %s',TrialData(i).Info.Condition,TrialData(i).Info.Trial);
        else
            titlename = sprintf('%s',TrialData(i).Info.Trial);
        end
        title(titlename);
    end
end