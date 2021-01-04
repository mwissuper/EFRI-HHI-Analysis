% Plot time series data for an individual subject. Plot only significant
% predictors from regression. Tension is positive for sensor in all
% directions. However, need to negate force in all directions when plotting
% force on POB.

clear; clc; %close all;

subj_array_force = [3:5 8:13];
s2 = [5 1 11 8.5]; % paper size
plotind = 0;
numcols = 2; numrows = 4;
colors(1,:) = [0.00,0.45,0.74]; % nice blue
colors(2,:) = [0.85,0.33,0.10]; % nice red
colors(3,:) = [0.47,0.67,0.19]; % nice green

for subj = 8%subj_array_force
    filename = sprintf('HHI2017_%i.mat',subj);
    load(filename);

    % Cycle through all trials and plot time series data for actual force,
    % predicted force, marker data, etc.

    %% ML/x dir L col

    for i = 1:length(TrialData)
       
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') %strcmp(TrialData(i).Info.Condition,'Solo Ground') %|| strcmp(TrialData(i).Info.Condition,'Assist Ground') 
%             figure;
            
            %% ML dir
                    
            plotind = 1; 
            
            % Plot POB Torso disp (sway)
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind),hold on;
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.torso(3:end,1)-TrialData(i).Results.beamMidline,'color',colors(3,:));
            % Labels and formatting
            ylabel('xb (m)'); box off; set(gca,'tickdir','out')
            
            % Plot POB Torso vel
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind),hold on;
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.vTorso(2:end,1),'color',colors(2,:))
            % Labels and formatting
            ylabel('xbdot (m/s)'); box off; set(gca,'tickdir','out');
            
            % Plot IP vel
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind),hold on;
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.vTorso(2:end,1),'color',colors(2,:))
            % Labels and formatting
            ylabel('xbdot (m/s)'); box off; set(gca,'tickdir','out');

%             % Plot power at interaction point (F x vFIN). 
%             subplot(numrows,numcols,plotind),hold on;
%             plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntPower(:,1));
%             hline(0,'k');
%             % Labels and formatting
%             ylabel('Int. Pt. Power ML (W)'); box off; set(gca,'tickdir','out');
%             titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Condition,TrialData(i).Info.Trial);
%             title(titlename);
%             
%             % Calculate scaling factors for acc and vel based on CLAV
%             temp = TrialData(i).Results.CLAV(3:end,1)-TrialData(i).Results.CLAV(1,1);
%             SFax = max(abs(temp))/max(abs(TrialData(i).Results.aCLAV(:,1)));
%             SFvx = max(abs(temp))/max(abs(TrialData(i).Results.vCLAV(2:end,1)));
            
            % Plot actual force
            
            %, modeled force, and components of clav model
            % scaled by regression coeff's. Must negate sign on Force
            % values since fitted to Force on POB
%             plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind)
%             plotind = plotind + numcols; 
            clear temp;
            temp.c = TrialData(i).Results.cx_torso;
            temp.c(isnan(temp.c)) = 0;
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(3:end,1)-temp.c(4),'k-','color',0.75.*[1 1 1],'linewidth',2),hold on;
%             m = [TrialData(i).Results.aTorso(:,1) TrialData(i).Results.vTorso(2:end,1) TrialData(i).Results.torso(3:end,1)-TrialData(i).Results.beamMidline];
%             plot(TrialData(i).Results.time(3:end),m*temp.c(1:3),'k');
%             % Plot modeled force components 
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cx_torso(1).*TrialData(i).Results.aTorso(:,1),'color',colors(1,:));
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cx_torso(2).*TrialData(i).Results.vTorso(2:end,1),'color',colors(2,:));
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cx_torso(3).*(TrialData(i).Results.torso(3:end,1)-TrialData(i).Results.beamMidline),'color',colors(3,:));
            % Labels and formatting
            ylabel('Force (N)'); box off; set(gca,'tickdir','out');
%             titlename = sprintf('R^2=%.2f m:%.2f b:%.2f k:%.2f',TrialData(i).Results.rsqx_torso,TrialData(i).Results.cx_torso(1),TrialData(i).Results.cx_torso(2),TrialData(i).Results.cx_torso(3));
%             legend('F','Fmodel','Fmass','Fdamper','Fspring','orientation','horizontal','location','northoutside');
            xlabel('Time (s)');
%             title(titlename);
            
%             % Plot POB Torso acc
%             plotind = plotind + numcols; 
%             subplot(numrows,numcols,plotind),hold on;
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.aTorso(:,1),'color',colors(1,:))
%             % Labels and formatting
%             ylabel('xbddot (m/s^2)'); box off; set(gca,'tickdir','out');
           
                   
            %% Vert dir
            
            plotind = 2;
            % Plot power at interaction point (F x vFIN) 
            subplot(numrows,numcols,plotind),hold on;
            plot(TrialData(i).Results.time(2:end),TrialData(i).Results.IntPower(:,3));
            hline(0,'k');
            % Labels and formatting
            ylabel('Int. Pt. Power Vert. (W)'); box off; set(gca,'tickdir','out');
            
            % Calculate scaling factors for acc and vel based on CLAV
            temp = TrialData(i).Results.CLAV(3:end,1)-TrialData(i).Results.CLAV(1,1);
            SFaz = max(abs(temp))/max(abs(TrialData(i).Results.aCLAV(:,3)));
            SFvz = max(abs(temp))/max(abs(TrialData(i).Results.vCLAV(2:end,3)));
                                   
            % Plot actual force, modeled force, and components of model
            % scaled by regression coeff's for clav
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind)
            clear temp;
            temp.c = TrialData(i).Results.cz_torso;
            temp.c(isnan(temp.c)) = 0;
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(3:end,3),'k-','color',0.75.*[1 1 1],'linewidth',2),hold on;
            m = [TrialData(i).Results.aCLAV(:,3) TrialData(i).Results.vCLAV(2:end,3) TrialData(i).Results.CLAV(3:end,3)-TrialData(i).Results.CLAV(1,3) ones(size(TrialData(i).Results.aCLAV(:,3)))];
            plot(TrialData(i).Results.time(3:end),m*temp.c,'k');
            % Plot modeled force components 
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cz_torso(1).*TrialData(i).Results.aCLAV(:,3),'color',colors(1,:));
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cz_torso(2).*TrialData(i).Results.vCLAV(2:end,3),'color',colors(2,:));
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cz_torso(3).*(TrialData(i).Results.CLAV(3:end,3)-TrialData(i).Results.CLAV(1,3)),'color',colors(3,:));
            % Labels and formatting
            ylabel('F Vert (N)'); box off; set(gca,'tickdir','out');
            titlename = sprintf('R^2=%.2f m:%.2f b:%.2f k:%.2f',TrialData(i).Results.rsqz_torso,TrialData(i).Results.cz_torso(1),TrialData(i).Results.cz_torso(2),TrialData(i).Results.cz_torso(3));
            legend('F','Fmodel','F acc clav','F vel clav','F disp clav','orientation','horizontal','location','northoutside');
            title(titlename);
            
            % Plot POB Clavicle state
            plotind = plotind + numcols; 
            subplot(numrows,numcols,plotind),hold on;
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.aCLAV(:,3).*SFaz,'color',colors(1,:))
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.vCLAV(2:end,3).*SFvz,'color',colors(2,:))
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.CLAV(3:end,3)-nanmean(TrialData(i).Results.CLAV(3:end,3)),'color',colors(3,:));
            % Labels and formatting
            ylabel('Clav Vert'); box off; set(gca,'tickdir','out');
            legend('acc','vel','disp','orientation','horizontal','location','northoutside');
            
            %% Save one pdf per subject
            set(gcf,'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
            pdfname = sprintf('HHI%i_modelF_power_SoloGround',subj);
            set(gcf, 'Color', 'w'); 
            export_fig(gcf,pdfname,'-pdf','-append')
            close all;
   
        end
    end
end

%% Cycle through all trials and plot time series data for force and vel (check work)

plotind = 0;
plotrow1 = 0;
plotrow2 = 0;
plotrow3 = 0;
plotrow4 = 0;
for i = 1:length(TrialData)
    if strcmp(TrialData(i).Info.Condition,'Solo Beam') || strcmp(TrialData(i).Info.Condition,'Assist Beam') %|| strcmp(TrialData(i).Info.Condition,'Solo Ground') || strcmp(TrialData(i).Info.Condition,'Assist Ground')
        plotind = plotind + 1;
        if strcmp(TrialData(i).Info.Condition,'Solo Beam')  
            titlename = sprintf('T%i, Solo Beam, rho = %.2f',i,r);
        elseif strcmp(TrialData(i).Info.Condition,'Assist Beam')
            titlename = sprintf('T%i, Asst Beam, rho = %.2f',i,r);
        end
        if plotind == 1
            titlename = sprintf('HHI %i %s',subj,titlename);
        end
        if plotind >= 1
            subplot(2,numcols,plotind),[ax,h1,h2] = plotyy(TrialData(i).Results.workInt,TrialData(i).Results.legObliq,TrialData(i).Results.time,TrialData(i).Results.thoraxObliq); 
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
