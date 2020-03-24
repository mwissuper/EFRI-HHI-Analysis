% Plot time series data for an individual subject to check against trial and subj
% means
clear; clc; close all;

analyze_force = 1;

% Load the concatenated group stats for late trials
if analyze_force == 1
    subj_array = [3:5 8:13];
    load('HHI2017_LateStats_force_MW.mat'); 
    conds = {'Assist Ground','Assist Beam'};
else
    subj_array = 3:14;
    load('HHI2017_LateStats_MW.mat'); % More participants than force data
    conds = {'Solo Beam','Assist Beam'};
end

% Constants
s2 = [5 1 11 8.5]; % paper size
numcols = 2; % 1 for ML dir, 2 for Vert dir
indSubj = 0;

for subj = 3%subj_array
    subj = subj
    indSubj = indSubj + 1;
    filename = sprintf('HHI2017_%i.mat',subj); % Load the individual subject's time series data
    load(filename);
    
    if analyze_force == 1
        numrows = 5; % Number of metrics
    else
        numrows = 2; % 2 performance metrics based on time series data
    end

    for i = 1:length(TrialData)
       
        if strcmp(TrialData(i).Info.Condition,conds{1}) || strcmp(TrialData(i).Info.Condition,conds{2}) 
            figure;          
            
            if strcmp(TrialData(i).Info.Condition,conds{1})
                indCond = 1;
            else
                indCond = 2;
            end
            
            %% ML dir
            
            % Abs force
            plotind = 1; 
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time,abs(TrialData(i).Results.Forces(:,1)));
            % Labels and formatting
            ylabel('Abs F ML (N)'); box off; set(gca,'tickdir','out');
            % Get mean corresponding to current trial in LateGroup data
            trialNum = str2num(TrialData(i).Info.Trial(end-1:end));
            temp = LateGroup.Subject==subj & LateGroup.TrialNumber==trialNum;
            row = find(temp==1,1,'first');
            if ~isempty(row) % empty if early trial
                hline(LateGroup.meanFx(row),'k--');
            end
            
            titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
            title(titlename);
                    
            % Pos/neg force
            plotind = plotind + numcols;
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time,TrialData(i).Results.Forces(:,1));
            if ~isempty(row) % empty if early trial
                hline(LateGroup.meanPosFx(row),'g-');
                hline(LateGroup.meanNegFx(row),'r-');
            end
            hline(0,'k--');
            % Labels and formatting
            ylabel('F ML (N)'); box off; set(gca,'tickdir','out');
                                  
            % Vel of Int. Pt. (used for power)
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time(2:end),TrialData(i).Results.vCLAV(:,1));
            hline(0,'k--');
            % Labels and formatting
            ylabel('Int. Pt. vel ML (m/s)'); box off; set(gca,'tickdir','out');
                     
            % Pos/neg power
            plotind = plotind + numcols;
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntPower(:,1));
            if ~isempty(row) % empty if early trial
                hline(LateGroup.meanPosPowerIntPtX(row),'g-');
                hline(LateGroup.meanNegPowerIntPtX(row),'r-');
            end
            hline(0,'k--');
            % Labels and formatting
            ylabel('Power ML (J/s)'); box off; set(gca,'tickdir','out');
            
            % Abs power
            plotind = plotind + numcols;
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time(2:end),abs(TrialData(i).Results.IntPower(:,1)));
            if ~isempty(row) % empty if early trial
                hline(LateGroup.meanAbsPowerIntPtX(row),'k-');
            end
            hline(0,'k--');
            % Labels and formatting
            ylabel('Abs power ML (J/s)'); box off; set(gca,'tickdir','out');
            
            % Plot xcorr clavicle disp. to force
            temp.F = TrialData(i).Results.Forces(1,3:end)-TrialData(i).Results.cx(4);
            temp.disp = TrialData(i).Results.vCLAV(1,3:end));
            temp.lag = TrialData(i).Results.lagX;
            if temp.lag < 0
                timeShift = TrialData(i).Results.time(3:end+temp.lag);
                swayShift = temp.disp(-temp.lag+1:end);
                temp.FCorr = temp.F(1:end+temp.lag);
            elseif temp.lag > 0
                timeShift = TrialData(i).Results.time(3+temp.lag:end);
                swayShift = temp.disp(1:end-temp.lag);
                temp.FCorr = temp.F(1+temp.lag:end);
            else
                timeShift = TrialData(i).Results.time(3:end);
                swayShift = temp.disp;
                temp.FCorr = temp.F;
            end
            temp.p = polyfit(swayShift,temp.FCorr,1);            
            % Plot force and sway xcorr
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind),hold on;
            plot(TrialData(i).Results.time(3:end),temp.F,'k-','color',0.75.*[1 1 1],'linewidth',2),hold on;
            plot(timeShift,polyval(temp.p,swayShift),'color',colors(3,:));
            % Labels and formatting
            box off; set(gca,'tickdir','out');
            titlename = sprintf('r = %.2f, lag = %.2f',TrialData(i).Results.xcorrX,temp.lag/TrialData(i).Markers.samplerate);
            title(titlename);
            legend('F ML','Clav ML disp','orientation','horizontal','location','northoutside');
            xlabel('Time (s)');
            
            temp = [];
            
%             % Plot pos/neg work
%             plotind = plotind + numcols;
%             subplot(numrows,numcols,plotind)
%             plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntCumWork(:,1));
% %             hline(gd.meanPosWorkIntPtX(indSubj,indCond),'g-');
% %             hline(gd.meanNegWorkIntPtX(indSubj,indCond),'r-');
%             hline(0,'k--');
%             % Labels and formatting
%             ylabel('Cum. Work ML (J)'); box off; set(gca,'tickdir','out');
            
%             % Abs work
%             plotind = plotind + numcols; 
%             subplot(numrows,numcols,plotind)
%             plot(TrialData(i).Results.time(2:end),abs(TrialData(i).Results.IntCumWork(:,1)));
%             hline(gd.meanWorkIntPtX(indSubj,indCond),'k-');
%             % Labels and formatting
%             ylabel('Abs Work ML (J)'); box off; set(gca,'tickdir','out');
            
            %% Vert dir
            
            % Abs force
            plotind = 2; 
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time,abs(TrialData(i).Results.Forces(:,3)));
            % Labels and formatting
            ylabel('Abs F Vert (N)'); box off; set(gca,'tickdir','out');
            % Get mean corresponding to current trial in LateGroup data
            trialNum = str2num(TrialData(i).Info.Trial(end-1:end));
            temp = LateGroup.Subject==subj & LateGroup.TrialNumber==trialNum;
            row = find(temp==1,1,'first');
            if ~isempty(row) % empty if early trial
                hline(LateGroup.meanFz(row),'k--');
            end
            
            % Plot pos/neg force
            plotind = plotind + numcols;
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time,TrialData(i).Results.Forces(:,3));
            if ~isempty(row) % empty if early trial
                hline(LateGroup.meanPosFz(row),'g-');
                hline(LateGroup.meanNegFz(row),'r-');
            end
            hline(0,'k--');
            % Labels and formatting
            ylabel('F Vert (N)'); box off; set(gca,'tickdir','out');
                        
            % Vel of Int. Pt. (used for power)
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time(2:end),TrialData(i).Results.vCLAV(:,3));
            hline(0,'k--');
            % Labels and formatting
            ylabel('Int. Pt. vel Vert (m/s)'); box off; set(gca,'tickdir','out');
            
            % Plot pos/neg power
            plotind = plotind + numcols;
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntPower(:,3));
            if ~isempty(row) % empty if early trial
                hline(LateGroup.meanPosPowerIntPtZ(row),'g-');
                hline(LateGroup.meanNegPowerIntPtZ(row),'r-');
            end
            hline(0,'k--');
            % Labels and formatting
            ylabel('Power Vert (J/s)'); box off; set(gca,'tickdir','out');
            
            % Abs power
            plotind = plotind + numcols;
            subplot(numrows,numcols,plotind)
            plot(TrialData(i).Results.time(2:end),abs(TrialData(i).Results.IntPower(:,3)));
            if ~isempty(row) % empty if early trial
                hline(LateGroup.meanAbsPowerIntPtZ(row),'k-');
            end
            hline(0,'k--');
            % Labels and formatting
            ylabel('Abs power Vert (J/s)'); box off; set(gca,'tickdir','out');
            
%             % Plot pos/neg work
%             plotind = plotind + numcols;
%             subplot(numrows,numcols,plotind)
%             plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntCumWork(:,3));
% %             hline(gd.meanPosWorkIntPtZ(indSubj,indCond),'g-');
% %             hline(gd.meanNegWorkIntPtZ(indSubj,indCond),'r-');
%             hline(0,'k--');
%             % Labels and formatting
%             ylabel('Cum. Work Vert (J)'); box off; set(gca,'tickdir','out');
            
%             % Abs work
%             plotind = plotind + numcols; 
%             subplot(numrows,numcols,plotind)
%             plot(TrialData(i).Results.time(2:end),abs(TrialData(i).Results.IntCumWork(:,3)));
%             hline(gd.meanWorkIntPtZ(indSubj,indCond),'k-');
%             % Labels and formatting
%             ylabel('Abs Work Vert (J)'); box off; set(gca,'tickdir','out');
            
            %% Save one pdf per subject
            set(gcf,'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
            if analyze_force == 1
                pdfname = sprintf('HHI%i_kinetics_time',subj);
            else
                pdfname = sprintf('HHI%i_kinem_time',subj);
            end
            set(gcf, 'Color', 'w'); 
            export_fig(gcf,pdfname,'-pdf','-append')
            close all;
   
        end
    end
end

%% Cycle through all trials and plot time series data for Pelvic and Thorax obliquitiy
% gpData = load('HHI2017_AllSubjectsLateStats.mat');
plotind = 0;
plotrow1 = 0;
plotrow2 = 0;
plotrow3 = 0;
plotrow4 = 0;
for i = 1:length(TrialData)
    if strcmp(TrialData(i).Info.Condition,'Solo Beam') || strcmp(TrialData(i).Info.Condition,'Assist Beam') %|| strcmp(TrialData(i).Info.Condition,'Solo Ground') || strcmp(TrialData(i).Info.Condition,'Assist Ground')
        % This is already calculated in getPointsStatsHHI_MW, but easier to
        % just redo calculations here instead of loading all group data and
        % then pulling out entries
        % Calculate correlation two obliquity angles to see if in or out of
        % phase. Use rho value only if statistically sig.
%         [rho p] = corr(TrialData(i).Results.pelvicObliq,TrialData(i).Results.thoraxObliq');
        plotind = plotind + 1;
        [rho p] = corr(TrialData(i).Results.legObliq',TrialData(i).Results.thoraxObliq');
        if p < 0.05
            r = rho;
        else
            r = nan;
        end
        if strcmp(TrialData(i).Info.Condition,'Solo Beam')  
            titlename = sprintf('T%i, Solo Beam, rho = %.2f',i,r);
%             plotrow1 = plotrow1 + 1;
%             plotind = plotrow1;
        elseif strcmp(TrialData(i).Info.Condition,'Assist Beam')
            titlename = sprintf('T%i, Asst Beam, rho = %.2f',i,r);
%             plotrow2 = plotrow2 + 1;
%             plotind = plotrow2+numcols;
%         elseif strcmp(TrialData(i).Info.Condition,'Solo Ground')
%             titlename = sprintf('Solo Gd, rho = %.2f',r);
%             plotrow3 = plotrow3 + 1;
%             plotind = plotrow3+2*numcols;
%         elseif strcmp(TrialData(i).Info.Condition,'Assist Ground')
%             titlename = sprintf('Asst Gd, rho = %.2f',r);
%             plotrow4 = plotrow4 + 1;
%             plotind = plotrow4+3*numcols;
        end
        if plotind == 1
            titlename = sprintf('HHI %i %s',subj,titlename);
        end
        % Plot two obliquity angles auto-scaled to see if time correlation
%         subplot(4,numcols,plotind),[ax,h1,h2] = plotyy(TrialData(i).Results.time,TrialData(i).Results.pelvicObliq,TrialData(i).Results.time,TrialData(i).Results.thoraxObliq); hold on; 
        if plotind >= 1
            subplot(2,numcols,plotind),[ax,h1,h2] = plotyy(TrialData(i).Results.time,TrialData(i).Results.legObliq,TrialData(i).Results.time,TrialData(i).Results.thoraxObliq); 
            set(ax(1),'ycolor','b'); set(h1,'color','b'); set(ax(2),'ycolor','g'); set(h2,'color','g');
            title(titlename);
        end
        if plotind == 1 || plotind == numcols+1 || plotind == 2*numcols+1 || plotind == 3*numcols+1
%             yyaxis left; a = ylabel('Pelvic obliq. (deg)'); set(a,'color','b')
            yyaxis left; a = ylabel('Leg obliq. (deg)'); set(a,'color','b')
        elseif plotind == numcols || plotind == 2*numcols || plotind == 3*numcols || plotind == 4*numcols
            yyaxis right; a = ylabel('Thorax obliq. (deg)'); set(a,'color','g')
        end
    end
end
set(gcf,'Position',   1.0e+03 *[0.0037    0.9180    2.5513    0.420]);

%% Cycle through all trials and plot time series data for Clav sway and F vert
plotind = 0;
plotrow1 = 0;
plotrow2 = 0;
for i = 1:length(TrialData)
    if strcmp(TrialData(i).Info.Condition,'Assist Beam')  
        titlename = sprintf('Asst Beam, rho = %.2f',TrialData(i).Results.rFvertClav);
        plotrow1 = plotrow1 + 1;
        plotind = plotrow1;
    elseif strcmp(TrialData(i).Info.Condition,'Assist Ground')
        titlename = sprintf('Asst Gd, rho = %.2f',TrialData(i).Results.rFvertClav);
        plotrow2 = plotrow2 + 1;
        plotind = plotrow2+numcols;
    end
    if plotind == 1
        titlename = sprintf('HHI %i',subj);
    end
    if strcmp(TrialData(i).Info.Condition,'Assist Beam') || strcmp(TrialData(i).Info.Condition,'Assist Ground')
        subplot(2,numcols,plotind),[ax,h1,h2] = plotyy(TrialData(i).Results.time,TrialData(i).Results.beamerSway,TrialData(i).Results.time,TrialData(i).Results.Forces(:,3)); hold on; % Clav sway and vertical force
        set(ax(1),'ycolor','b');
        if plotind == 1 || plotind == numcols+1
            yyaxis left; a = ylabel('Clav sway (mm)'); set(a,'color','b')
        elseif plotind == numcols || plotind == 2*numcols
            yyaxis right; a = ylabel('Vert. Force (N)'); set(a,'color','g')
        end
        set(h1,'color','b');
        set(ax(2),'ycolor','g');
        set(h2,'color','g');
        title(titlename);
    end
end
set(gcf,'Position',   1.0e+03 *[0.0037    0.9180    2.5513    0.420]);
