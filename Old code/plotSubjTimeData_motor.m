% Plot time series data for an individual subject. Plot residual of force fit for
% periods of time where power is motor. Must negate signs when plotting
% force and power because looking at force on POB
clear; clc; close all;

subj_array_force = [3:5 8:13];
s2 = [5 1 8.5 11]; % paper size
plotind = 0;
numcols = 1; numrows = 3;
colors(1,:) = [0.00,0.45,0.74]; % nice blue
colors(2,:) = [0.85,0.33,0.10]; % nice red
colors(3,:) = [0.47,0.67,0.19]; % nice green

for subj = 4%subj_array_force
    filename = sprintf('HHI2017_%i.mat',subj);
%     filename = 'HHI2017_8_FbiasAsstSolo.mat';
    load(filename);

    % Cycle through all trials and plot time series data for actual force,
    % predicted force, marker data, etc.

    %% ML/x dir L col

    for i = 1:length(TrialData)
       
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') %|| strcmp(TrialData(i).Info.Condition,'Assist Ground') 
            figure;
            
            %% ML dir
            plotind = 1; 
            
            % Find where power is positive or negative (motor or brake)
            temp.indNeg = find(TrialData(i).Results.IntPower(3:end,1) < 0); % Only look at time indices where we have acc data. This index is relative to acc data vector
           
            % Plot power at interaction point (F x vFIN). Corrected signs
            % in mainWorkPowerAnalysisMW.m s.t. all positive power
            % corresponds to motor in all directions
            subplot(numrows,numcols,plotind),hold on;
            plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntPower(:,1),'k-');
            temp.posPower = TrialData(i).Results.IntPower(:,1);
            temp.posPower(temp.indNeg+2) = nan;
            plot(TrialData(i).Results.time(2:end),temp.posPower,'k-','linewidth',2);
            hline(0,'k');
            % Labels and formatting
            ylabel('Int. Pt. Power ML (W)'); box off; set(gca,'tickdir','out');
            titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Condition,TrialData(i).Info.Trial);
            title(titlename);
            
            % Calculate scaling factors for acc and vel based on CLAV
            temp.range = TrialData(i).Results.CLAV(3:end,1)-TrialData(i).Results.CLAV(1,1);
            SFax = max(abs(temp.range))/max(abs(TrialData(i).Results.aCLAV(:,1)));
            SFvx = max(abs(temp.range))/max(abs(TrialData(i).Results.vCLAV(2:end,1)));
                      
            % Plot actual force, modeled force, and components of clav model
            % scaled by regression coeff's. 
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind)
            temp.c = TrialData(i).Results.cx_clav;
            temp.c(isnan(temp.c)) = 0;
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(3:end,1),'k-','color',0.75.*[1 1 1],'linewidth',2),hold on;
            m = [TrialData(i).Results.aCLAV(:,1) TrialData(i).Results.vCLAV(2:end,1) TrialData(i).Results.CLAV(3:end,1)-TrialData(i).Results.beamMidline ones(size(TrialData(i).Results.aCLAV(:,1)))];
            plot(TrialData(i).Results.time(3:end),m*temp.c,'k');
            % Plot modeled force components 
            temp.Fa = TrialData(i).Results.cx_clav(1).*TrialData(i).Results.aCLAV(:,1);
            temp.Fv = TrialData(i).Results.cx_clav(2).*TrialData(i).Results.vCLAV(2:end,1);
            temp.Fd = TrialData(i).Results.cx_clav(3).*(TrialData(i).Results.CLAV(3:end,1)-TrialData(i).Results.beamMidline);
            plot(TrialData(i).Results.time(3:end),temp.Fa,'color',colors(1,:));
            plot(TrialData(i).Results.time(3:end),temp.Fv,'color',colors(2,:));
            plot(TrialData(i).Results.time(3:end),temp.Fd,'color',colors(3,:));
            % Put in legend before plot more components
            legend('F','Fmodel','Fmass','Fdamper','Fspring','orientation','horizontal','location','northoutside');
            % Highlight force components where power IP is positive (motor)
            temp.Fa(temp.indNeg) = nan;
            temp.Fv(temp.indNeg) = nan;
            temp.Fd(temp.indNeg) = nan;
            plot(TrialData(i).Results.time(3:end),temp.Fa,'color',colors(1,:),'linewidth',2);
            plot(TrialData(i).Results.time(3:end),temp.Fv,'color',colors(2,:),'linewidth',2);
            plot(TrialData(i).Results.time(3:end),temp.Fd,'color',colors(3,:),'linewidth',2);
            
            % Labels and formatting
            ylabel('Mediolateral Force (N)'); box off; set(gca,'tickdir','out');
            titlename = sprintf('R^2=%.2f m:%.2f b:%.2f k:%.2f',TrialData(i).Results.rsqx_clav,TrialData(i).Results.cx_clav(1),TrialData(i).Results.cx_clav(2),TrialData(i).Results.cx_clav(3));
            xlabel('Time (s)');
            title(titlename);
            
            % Plot residual of model fit and highlight sections where power
            % is positive. Show value of correlation residual to power.
            temp.resid = TrialData(i).Results.Forces(3:end,1) - m*temp.c; % Must negate sign to get Force on POB and then compare to model
            % Calc corr power to F resid
            [temp.rho, temp.p] = corr(TrialData(i).Results.IntPower(2:end,1),temp.resid);
            if temp.p < 0.05
                titlename = sprintf('Corr. power vs. Fresid = %.2f',temp.rho); 
            else
                titlename = sprintf('Corr. p = %.2f',temp.p); 
            end
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind),hold on;
            plot(TrialData(i).Results.time(3:end),temp.resid,'k-'),hold on;
            temp.resid(temp.indNeg) = nan;
            plot(TrialData(i).Results.time(3:end),temp.resid,'k-','linewidth',2)
            hline(0,'k--');
            % Labels and formatting
            ylabel('F (N)'); box off; set(gca,'tickdir','out');
            legend('model residual','model residual p > 0','orientation','horizontal','location','northoutside');
            temp = [];
            title(titlename);
            
%             % Plot clavicle and rFin velocity to check signs on power
%             plotind = plotind + numcols; 
%             subplot(numrows,numcols,plotind),hold on;
%             plot(TrialData(i).Results.time(2:end),TrialData(i).Results.vCLAV(:,1))
%             plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntPtVel(:,1)) 
%             % Labels and formatting
%             ylabel('Lateral velocity (m/s)'); box off; set(gca,'tickdir','out');
%             legend('clav','IP','orientation','horizontal','location','northoutside');
%             
%             
            %% Save one pdf per subject
            set(gcf,'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
            pdfname = sprintf('HHI%i_modelF_motor_AssistBeam',subj);
            set(gcf, 'Color', 'w'); 
            export_fig(gcf,pdfname,'-pdf','-append')
            close all;
   
        end
    end
end

