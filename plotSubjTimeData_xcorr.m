% Plot time series data for an individual subject to check against trial and subj
% means. ML direction only
clear; clc; close all;

% Load the concatenated group stats for late trials. Need force data
subj_array = [3:5 8:13];
load('HHI2017_LateStats_force_MW.mat'); 
conds = {'Assist Ground','Assist Beam'};

% Constants
s2 = [5 1 8.5 11]; % portrait paper size. later print it out in booklet style with two pages per paper page
numcols = 1; numrows = 4; 
indSubj = 0;

for subj = subj_array
    subj = subj
    indSubj = indSubj + 1;
    filename = sprintf('HHI2017_%i.mat',subj); % Load the individual subject's time series data
    load(filename);
    pdfnames{1} = sprintf('HHI%i_IP_Clav_time_Ground',subj);
    pdfnames{2} = sprintf('HHI%i_IP_Clav_time_Beam',subj);
    
    for j = 1:length(conds)       
        for i = 1:length(TrialData)
            if strcmp(TrialData(i).Info.Condition,conds{j})
                figure;          

                if strcmp(TrialData(i).Info.Condition,conds{1})
                    indCond = 1;
                else
                    indCond = 2;
                end

                %% F at IP overlaid with vClav, show max xcorr and lag
                plotind = 1; 
                subplot(numrows,numcols,plotind)
                yyaxis left
                plot(TrialData(i).Results.time(2:end),TrialData(i).Results.Forces(2:end,1)),ylabel('Force (N)');
                yyaxis right
                plot(TrialData(i).Results.time(2:end),TrialData(i).Results.vCLAV(:,1)),ylabel('v Clav (m/s)');
                % Labels and formatting
                box off; set(gca,'tickdir','out');
                titlename = sprintf('HHI%i %s %s, maxXcorr = %.2f, lag = %.2f',...
                    subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition,...
                    TrialData(i).Results.xcorrFIPvClavX,TrialData(i).Results.lagFIPvClavX);
                title(titlename);

                %% vIP overlaid with vClav, show max xcorr and lag
                plotind = plotind + 1; 
                subplot(numrows,numcols,plotind)
                plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntPtVel(:,1)),hold on
                plot(TrialData(i).Results.time(2:end),TrialData(i).Results.vCLAV(:,1))
                % Labels and formatting
                ylabel('vel (m/s)'),box off; set(gca,'tickdir','out');
                legend('IP','Clav');
                titlename = sprintf('maxXcorr = %.2f, lag = %.2f',...
                    TrialData(i).Results.xcorrvIPvClavX,TrialData(i).Results.lagvIPvClavX);
                title(titlename);

                %% arm length of POB
                plotind = plotind + 1; 
                subplot(numrows,numcols,plotind)
                plot(TrialData(i).Results.time,TrialData(i).Results.armPOBX)
                % Labels and formatting
                box off; set(gca,'tickdir','out');
                ylabel('ML dist. IP to Clav (m)');
                                
                %% Find where power is positive or negative (motor or brake)
                temp.indNeg = find(TrialData(i).Results.IntPower(3:end,1) < 0); % Only look at time indices where we have acc data. This index is relative to acc data vector

                % Plot power at interaction point (F x vFIN). Corrected signs
                % in mainWorkPowerAnalysisMW.m s.t. all positive power
                % corresponds to motor in all directions
                plotind = plotind + 1; 
                subplot(numrows,numcols,plotind),hold on;
                plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntPower(:,1),'k-');
                temp.posPower = TrialData(i).Results.IntPower(:,1);
                temp.posPower(temp.indNeg+2) = nan;
                plot(TrialData(i).Results.time(2:end),temp.posPower,'k-','linewidth',2);
                hline(0,'k');
                % Labels and formatting
                ylabel('IP Power ML (W)'); box off; set(gca,'tickdir','out');
                title(titlename);
                
                xlabel('Time (s)');

    %             % Plot xcorr clavicle disp. to force
    %             temp.F = TrialData(i).Results.Forces(1,3:end)-TrialData(i).Results.cx(4);
    %             temp.disp = TrialData(i).Results.vCLAV(1,3:end));
    %             temp.lag = TrialData(i).Results.lagX;
    %             if temp.lag < 0
    %                 timeShift = TrialData(i).Results.time(3:end+temp.lag);
    %                 swayShift = temp.disp(-temp.lag+1:end);
    %                 temp.FCorr = temp.F(1:end+temp.lag);
    %             elseif temp.lag > 0
    %                 timeShift = TrialData(i).Results.time(3+temp.lag:end);
    %                 swayShift = temp.disp(1:end-temp.lag);
    %                 temp.FCorr = temp.F(1+temp.lag:end);
    %             else
    %                 timeShift = TrialData(i).Results.time(3:end);
    %                 swayShift = temp.disp;
    %                 temp.FCorr = temp.F;
    %             end
    %             temp.p = polyfit(swayShift,temp.FCorr,1);            
    %             % Plot force and sway xcorr
    %             plotind = plotind + numcols; 
    %             subplot(numrows,numcols,plotind),hold on;
    %             plot(TrialData(i).Results.time(3:end),temp.F,'k-','color',0.75.*[1 1 1],'linewidth',2),hold on;
    %             plot(timeShift,polyval(temp.p,swayShift),'color',colors(3,:));
    %             % Labels and formatting
    %             box off; set(gca,'tickdir','out');
    %             titlename = sprintf('r = %.2f, lag = %.2f',TrialData(i).Results.xcorrX,temp.lag/TrialData(i).Markers.samplerate);
    %             title(titlename);
    %             legend('F ML','Clav ML disp','orientation','horizontal','location','northoutside');
    %             xlabel('Time (s)');
    %             
    %             temp = [];

                %% Save one pdf per subject
                set(gcf,'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
                pdfname = pdfnames{j};
                set(gcf, 'Color', 'w'); 
                export_fig(gcf,pdfname,'-pdf','-append')
                close all;

            end
        end
    end
end
