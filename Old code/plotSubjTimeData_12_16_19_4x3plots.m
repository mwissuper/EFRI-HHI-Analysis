% Plot time series data for an individual subject
clear; clc; close all;

subj_array_force = [3:5 8:13];
for subj = subj_array_force
    filename = sprintf('HHI2017_%i.mat',subj);
    load(filename);
    numcols = 3; numrows = 4;

    % Cycle through all trials and plot time series data for actual force, predicted force, each force component multiplied by regression coefficient

    %% x dir
    plotind1 = 0; plotind2 = 0; plotind3 = 0; plotind4 = 0;
    h1 = figure(1); h2 = figure(2); h3 = figure(3); h4 = figure(4);
    for i = 1:length(TrialData)
        % Plot forces
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') 
            figure(1);
            plotind1 = plotind1 + 1; 
            subplot(numrows,numcols,plotind1)
        elseif strcmp(TrialData(i).Info.Condition,'Assist Ground') 
            figure(2);
            plotind2 = plotind2 + 1;  
            subplot(numrows,numcols,plotind2)
        end
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') || strcmp(TrialData(i).Info.Condition,'Assist Ground') 
            % Plot (actual force - mean), modeled force
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(1,3:end)-TrialData(i).Results.cx(4),'k--'),hold on;
            m = [TrialData(i).Results.IntPtAcc(1,:)' TrialData(i).Results.IntPtVel(1,2:end)' TrialData(i).Results.IntPt(1,3:end)'-mean(TrialData(i).Results.IntPt(1,3:end))];
            plot(TrialData(i).Results.time(3:end),m*TrialData(i).Results.cx(1:3),'k');
            % Plot modeled force components 
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cx(1).*TrialData(i).Results.IntPtAcc(1,:));
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cx(2).*TrialData(i).Results.IntPtVel(1,2:end));
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cx(3).*(TrialData(i).Results.IntPt(1,3:end)-mean(TrialData(i).Results.IntPt(1,3:end))));
            % Labels and formatting
            xlabel('Time (s)'),ylabel('F ML (N)'); box off; set(gca,'tickdir','out');
            if (plotind1 == 2 && strcmp(TrialData(i).Info.Condition,'Assist Beam')) || (plotind2 == 2 &&  strcmp(TrialData(i).Info.Condition,'Assist Ground') )
                titlename = sprintf('HHI %i %s %s R^2=%.2f',subj,TrialData(i).Info.Condition,TrialData(i).Info.Trial,TrialData(i).Results.rsqx);
                legend('F','Fmodel','F acc','F vel','F disp','orientation','horizontal','location','south');
            else
                titlename = sprintf('%s R^2=%.2f',TrialData(i).Info.Trial,TrialData(i).Results.rsqx);
            end
            title(titlename);
        end

        % Plot force difference with clavicle marker state
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') || strcmp(TrialData(i).Info.Condition,'Assist Ground') 
            % Scaling factors for plotting clavicle state (rough method of take max of
            % F with mean removed and divide by max of clav state with mean removed)
            Fmrem = TrialData(i).Results.Forces(1,3:end)'-TrialData(i).Results.cx(4);
            SFa = max(Fmrem)./max(TrialData(i).Results.aCLAV(1,3:end)-nanmean(TrialData(i).Results.aCLAV(1,3:end)))';
            SFv = max(Fmrem)./max(TrialData(i).Results.vCLAV(1,3:end)-nanmean(TrialData(i).Results.vCLAV(1,3:end)))';
            SFd = max(Fmrem)./max(TrialData(i).Results.CLAV(1,3:end)-nanmean(TrialData(i).Results.CLAV(1,3:end)))';
        end
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') 
            figure(3);
            plotind3 = plotind3 + 1;  
            subplot(numrows,numcols,plotind3)
        elseif strcmp(TrialData(i).Info.Condition,'Assist Ground') 
            figure(4);
            plotind4 = plotind4 + 1;  
            subplot(numrows,numcols,plotind4)
        end
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') || strcmp(TrialData(i).Info.Condition,'Assist Ground') 
            % Plot (actual force - mean), modeled force, and difference in
            % black
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(1,3:end)-TrialData(i).Results.cx(4),'k--'),hold on;
%             m = [TrialData(i).Results.IntPtAcc(1,:)' TrialData(i).Results.IntPtVel(1,2:end)' TrialData(i).Results.IntPt(1,3:end)'-mean(TrialData(i).Results.IntPt(1,3:end))];
%             plot(TrialData(i).Results.time(3:end),m*TrialData(i).Results.cx(1:3),'k');
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(1,3:end)'-TrialData(i).Results.cx(4)-m*TrialData(i).Results.cx(1:3),'k-'); hold on;
            % Plot POB clavicle state multiplied by SF to fit on yaxes
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.aCLAV(1,:)'.*SFa);
            plot(TrialData(i).Results.time(3:end),TrialData(i).Results.vCLAV(1,2:end)'.*SFv);
            plot(TrialData(i).Results.time(3:end),(TrialData(i).Results.CLAV(1,3:end)'-nanmean(TrialData(i).Results.CLAV(1,3:end))).*SFd);
            % Labels and formatting
            xlabel('Time (s)'),ylabel('F ML (N)'); box off; set(gca,'tickdir','out');
            if (plotind3 == 2 && strcmp(TrialData(i).Info.Condition,'Assist Beam')) || (plotind4 == 2 && strcmp(TrialData(i).Info.Condition,'Assist Ground') )
                titlename = sprintf('HHI %i %s %s',subj,TrialData(i).Info.Condition,TrialData(i).Info.Trial);
%                 legend('F','Fmodel','F - Fmodel','Clav acc','Clav vel','Clav disp','orientation','horizontal','location','south');
                legend('F','F - Fmodel','Clav acc','Clav vel','Clav disp','orientation','horizontal','location','south');
            else
                titlename = sprintf('%s',TrialData(i).Info.Trial);
            end
            title(titlename);
        end
    end

    % Save fig files and also print to pdf

    s2 = [5 1 11 8.5];

    set(figure(1),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
    figname1 = sprintf('HHI%i_F_Fmodel_compX_AB.fig',subj);
    saveas(h1,figname1,'fig');
    pdfname1 = figname1(1:end-4);
    print(pdfname1,'-dpdf','-fillpage')

    set(figure(2),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
    figname2 = sprintf('HHI%i_F_Fmodel_compX_AG.fig',subj);
    saveas(h2,figname2,'fig');
    pdfname2 = figname2(1:end-4);
    print(pdfname2,'-dpdf','-fillpage')

    set(figure(3),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
    figname3 = sprintf('HHI%i_F_Clav_compX_AB.fig',subj);
    saveas(h3,figname3,'fig');
    pdfname3 = figname3(1:end-4);
    print(pdfname3,'-dpdf','-fillpage')

    set(figure(4),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
    figname4 = sprintf('HHI%i_F_Clav_compX_AG.fig',subj);
    saveas(h4,figname4,'fig');
    pdfname4 = figname4(1:end-4);
    print(pdfname4,'-dpdf','-fillpage')
    
    close all

%     %% z dir
%     plotind1 = 0; plotind2 = 0; plotind3 = 0; plotind4 = 0;
%     h1 = figure(1); h2 = figure(2); h3 = figure(3); h4 = figure(4);
%     for i = 1:length(TrialData)
%         % Plot forces
%         if strcmp(TrialData(i).Info.Condition,'Assist Beam') 
%             figure(1);
%             plotind1 = plotind1 + 1;  
%             subplot(numrows,numcols,plotind1)
%         elseif strcmp(TrialData(i).Info.Condition,'Assist Ground') 
%             figure(2);
%             plotind2 = plotind2 + 1;  
%             subplot(numrows,numcols,plotind2)
%         end
%         if strcmp(TrialData(i).Info.Condition,'Assist Beam') || strcmp(TrialData(i).Info.Condition,'Assist Ground') 
%             % Plot (actual force - mean), modeled force
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(3,3:end)-TrialData(i).Results.cz(4),'k--'),hold on;
%             m = [TrialData(i).Results.IntPtAcc(3,:)' TrialData(i).Results.IntPtVel(3,2:end)' TrialData(i).Results.IntPt(3,3:end)'-mean(TrialData(i).Results.IntPt(3,3:end))];
%             plot(TrialData(i).Results.time(3:end),m*TrialData(i).Results.cz(1:3),'k');
%             % Plot modeled force components 
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cz(1).*TrialData(i).Results.IntPtAcc(3,:));
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cz(2).*TrialData(i).Results.IntPtVel(3,2:end));
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.cz(3).*(TrialData(i).Results.IntPt(3,3:end)-mean(TrialData(i).Results.IntPt(3,3:end))));
%             % Labels and formatting
%             xlabel('Time (s)'),ylabel('F Vert (N)'); box off; set(gca,'tickdir','out');
%             if (plotind1 == 2 && strcmp(TrialData(i).Info.Condition,'Assist Beam')) || (plotind2 == 2 &&  strcmp(TrialData(i).Info.Condition,'Assist Ground') )
%                 titlename = sprintf('HHI %i %s %s R^2=%.2f',subj,TrialData(i).Info.Condition,TrialData(i).Info.Trial,TrialData(i).Results.rsqz);
%                 legend('F','Fmodel','F acc','F vel','F disp','orientation','horizontal','location','south');
%             else
%                 titlename = sprintf('%s R^2=%.2f',TrialData(i).Info.Trial,TrialData(i).Results.rsqz);
%             end
%             title(titlename);
%         end
% 
%         % Plot forces with clavicle marker state
%         if strcmp(TrialData(i).Info.Condition,'Assist Beam') || strcmp(TrialData(i).Info.Condition,'Assist Ground') 
%             % Scaling factors for plotting clavicle state (rough method of take max of
%             % F with mean removed and divide by max of clav state with mean removed)
%             Fmrem = TrialData(i).Results.Forces(3,3:end)'-TrialData(i).Results.cz(4);
%             SFa = max(Fmrem)./max(TrialData(i).Results.aCLAV(3,3:end)-nanmean(TrialData(i).Results.aCLAV(3,3:end)))';
%             SFv = max(Fmrem)./max(TrialData(i).Results.vCLAV(3,3:end)-nanmean(TrialData(i).Results.vCLAV(3,3:end)))';
%             SFd = max(Fmrem)./max(TrialData(i).Results.CLAV(3,3:end)-nanmean(TrialData(i).Results.CLAV(3,3:end)))';
%         end
%         if strcmp(TrialData(i).Info.Condition,'Assist Beam') 
%             figure(3);
%             plotind3 = plotind3 + 1;  
%             subplot(numrows,numcols,plotind3)
%         elseif strcmp(TrialData(i).Info.Condition,'Assist Ground') 
%             figure(4);
%             plotind4 = plotind4 + 1;  
%             subplot(numrows,numcols,plotind4)
%         end
%         if strcmp(TrialData(i).Info.Condition,'Assist Beam') || strcmp(TrialData(i).Info.Condition,'Assist Ground') 
% %             % Plot (actual force - mean), modeled force, and difference in
% %             % black
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(3,3:end)-TrialData(i).Results.cz(4),'k--'),hold on;
% %             m = [TrialData(i).Results.IntPtAcc(3,:)' TrialData(i).Results.IntPtVel(3,2:end)' TrialData(i).Results.IntPt(3,3:end)'-mean(TrialData(i).Results.IntPt(3,3:end))];
% %             plot(TrialData(i).Results.time(3:end),m*TrialData(i).Results.cz(1:3),'k');
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.Forces(3,3:end)'-TrialData(i).Results.cz(4)-m*TrialData(i).Results.cz(1:3),'k-'); hold on;
%             % Plot POB clavicle state multiplied by SF to fit on yaxes
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.aCLAV(3,:)'.*SFa);
%             plot(TrialData(i).Results.time(3:end),TrialData(i).Results.vCLAV(3,2:end)'.*SFv);
%             plot(TrialData(i).Results.time(3:end),(TrialData(i).Results.CLAV(3,3:end)'-nanmean(TrialData(i).Results.CLAV(3,3:end))).*SFd);
%             % Labels and formatting
%             xlabel('Time (s)'),ylabel('F Vert (N)'); box off; set(gca,'tickdir','out');
%             if (plotind3 == 2 && strcmp(TrialData(i).Info.Condition,'Assist Beam')) || (plotind4 == 2 && strcmp(TrialData(i).Info.Condition,'Assist Ground') )
%                 titlename = sprintf('HHI %i %s %s',subj,TrialData(i).Info.Condition,TrialData(i).Info.Trial);
%                 legend('F','F - Fmodel','Clav acc','Clav vel','Clav disp','orientation','horizontal','location','south');
%             else
%                 titlename = sprintf('%s',TrialData(i).Info.Trial);
%             end
%             title(titlename);
%         end
%     end
% 
%     % Save fig files and also print to pdf
% 
%     s2 = [5 1 11 8.5];
% 
%     set(figure(1),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
%     figname1 = sprintf('HHI%i_F_Fmodel_compZ_AB.fig',subj);
%     saveas(h1,figname1,'fig');
%     pdfname1 = figname1(1:end-4);
%     print(pdfname1,'-dpdf','-fillpage')
% 
%     set(figure(2),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
%     figname2 = sprintf('HHI%i_F_Fmodel_compZ_AG.fig',subj);
%     saveas(h2,figname2,'fig');
%     pdfname2 = figname2(1:end-4);
%     print(pdfname2,'-dpdf','-fillpage')
% 
%     set(figure(3),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
%     figname3 = sprintf('HHI%i_F_Clav_compZ_AB.fig',subj);
%     saveas(h3,figname3,'fig');
%     pdfname3 = figname3(1:end-4);
%     print(pdfname3,'-dpdf','-fillpage')
% 
%     set(figure(4),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
%     figname4 = sprintf('HHI%i_F_Clav_compZ_AG.fig',subj);
%     saveas(h4,figname4,'fig');
%     pdfname4 = figname4(1:end-4);
%     print(pdfname4,'-dpdf','-fillpage')

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
        subplot(2,numcols,plotind),[ax,h1,h2] = plotyy(TrialData(i).Results.time,TrialData(i).Results.beamerSway,TrialData(i).Results.time,TrialData(i).Results.Forces(3,:)); hold on; % Clav sway and vertical force
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
