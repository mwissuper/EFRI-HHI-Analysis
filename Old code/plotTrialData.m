% Plot individual trial data for a participant to understand force, power,
% and work

close all; clear all; clc;

subj = 3;

fname = sprintf('HHI2017_%i.mat',subj);
data = load(fname);

plotoption = 2; % 1 = plot work vs. time, 2 = plot force vs. time, 3 = plot vel RFIN vs. time
dir = 1; % 1 = ML or x direction, 2 = AP or y dir, 3 = vert or z
conds = {'Assist Ground','Assist Beam'};

% Plots of a given metric for all trials for one subj
plotind = 0;
numrows = 4; numcols = 5;
for i = 1:length(data.TrialData)
    if ~isempty(find(strcmp(data.TrialData(i).Info.Condition,conds)==1,1,'first'))
        plotind = plotind + 1;
        subplot(numrows,numcols,plotind)
        if plotoption == 1
            plot(data.TrialData(i).Results.time(2:end),data.TrialData(i).Results.IntCumWork(dir,:))
            xlabel('Time (s)'), ylabel('Work (J)'); hline(0,'k--');
            if dir == 1
                ylim([-1.1 .3]); 
            elseif dir == 3
                ylim([-.4 1]); 
            end
        elseif plotoption == 2
            plot(data.TrialData(i).Results.time,data.TrialData(i).Results.Forces(dir,:))
            xlabel('Time (s)'), ylabel('Force (N)'); hline(0,'k--');
            if dir == 1
                ylim([-4 17]); 
%             elseif dir == 3
%                 ylim([-.4 1]); 
            end
        end
        temp = data.TrialData(i).Info.Trial;
        trial = str2num(temp(end-1:end));
        if plotind == 1
            titlename = sprintf('S%i, Tr. %i, %s',subj,trial,data.TrialData(i).Info.Condition);
        else
            titlename = sprintf('Tr. %i, %s',trial,data.TrialData(i).Info.Condition);
        end
        title(titlename);
    end
end

%% Plots of velocity, force, and power at interaction point for a single trial
%% Lateral dir
subplot(4,1,1),plot(time(2:end),vFin(1,:)),xlabel('Time (s)'),ylabel('Lateral vFin (m/s)')
hline(0,'k--');
subplot(4,1,2),plot(time(2:end),Force(1,2:end)),xlabel('Time (s)'),ylabel('Lateral F (N)')
subplot(4,1,3),plot(time(2:end),intPt_power(1,:)),xlabel('Time (s)'),ylabel('Lateral Power (W)'),hold on;
a = find(intPt_power(1,:)>0);
b = find(intPt_power(1,:)<0);
hline(0,'k--');
hline(mean(intPt_power(1,a)),'r--');hline(mean(intPt_power(1,b)),'b--');
subplot(4,1,4),plot(time(2:end),intPt_cumWork(1,:)),xlabel('Time (s)'),ylabel('Lateral Cum Work (J)')
hline(0,'k--');
set(gcf,'outerposition',s);
%% Sagittal dir
subplot(3,1,1),plot(time(2:end),vFin(2,:)),xlabel('Time (s)'),ylabel('Sagittal vFin (m/s)')
subplot(3,1,2),plot(time(2:end),Force(2,2:end)),xlabel('Time (s)'),ylabel('Sagittal F (N)')
subplot(3,1,3),plot(time(2:end),intPt_power(2,:)),xlabel('Time (s)'),ylabel('Sagittal Power (W)')
hline(0,'k--');
set(gcf,'outerposition',s);
%% Vertical dir
subplot(4,1,1),plot(time(2:end),vFin(3,:)),xlabel('Time (s)'),ylabel('Vertical vFin (m/s)')
hline(0,'k--');
subplot(4,1,2),plot(time(2:end),Force(3,2:end)),xlabel('Time (s)'),ylabel('Vertical F (N)')
subplot(4,1,3),plot(time(2:end),intPt_power(3,:)),xlabel('Time (s)'),ylabel('Vertical Power (W)')
a = find(intPt_power(3,:)>0);
b = find(intPt_power(3,:)<0);
hline(0,'k--');
hline(mean(intPt_power(3,a)),'r--');hline(mean(intPt_power(3,b)),'b--');
subplot(4,1,4),plot(time(2:end),intPt_cumWork(3,:)),xlabel('Time (s)'),ylabel('Vertical Cum Work (W)')
hline(0,'k--');
set(gcf,'outerposition',s);

%% Plots of power and cum work at interaction point compared to markers on POB's body to see effect
%% Lateral dir
subplot(3,1,1),plot(time(2:end),vFin(1,:)),xlabel('Time (s)'),ylabel('Lateral vFin (m/s)')
hline(0,'k--');
subplot(3,1,2),plot(time(2:end),Force(1,2:end)),xlabel('Time (s)'),ylabel('Lateral F (N)')
subplot(3,1,3),plot(time(2:end),intPt_power(1,:)),xlabel('Time (s)'),ylabel('Lateral Power (W)')
set(gcf,'outerposition',s);
%% Sagittal dir
subplot(3,1,1),plot(time(2:end),vFin(2,:)),xlabel('Time (s)'),ylabel('Sagittal vFin (m/s)')
subplot(3,1,2),plot(time(2:end),Force(2,2:end)),xlabel('Time (s)'),ylabel('Sagittal F (N)')
subplot(3,1,3),plot(time(2:end),intPt_power(2,:)),xlabel('Time (s)'),ylabel('Sagittal Power (W)')
hline(0,'k--');
set(gcf,'outerposition',s);
%% Vertical dir
subplot(3,1,1),plot(time(2:end),vFin(3,:)),xlabel('Time (s)'),ylabel('Vertical vFin (m/s)')
hline(0,'k--');
subplot(3,1,2),plot(time(2:end),Force(3,2:end)),xlabel('Time (s)'),ylabel('Vertical F (N)')
subplot(3,1,3),plot(time(2:end),intPt_power(3,:)),xlabel('Time (s)'),ylabel('Vertical Power (W)')
hline(0,'k--');
set(gcf,'outerposition',s);

%%
title('Assist Ground');
%%
title('Assist Beam');
%% line
hline(0,'k--');