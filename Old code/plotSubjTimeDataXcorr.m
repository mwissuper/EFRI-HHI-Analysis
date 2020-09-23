% Plot time series data for an individual subject
clear; clc; close all;

subj_array_force = [3:5 8:13];
s2 = [5 1 11 8.5]; % paper size

for subj = subj_array_force
    filename = sprintf('HHI2017_%i.mat',subj);
    load(filename);
    fs = TrialData(1).Markers.samplerate;
    numcols = 2; numrows = 2;

    % Cycle through all trials and plot time series data for actual force
    % and sway in ML and Vert. dir's and combo's for xcorr range (-1s to
    % +1s)

    %% ML/x force L col, Vertical/z force R col

    for i = 1:length(TrialData)
       
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') || strcmp(TrialData(i).Info.Condition,'Assist Ground') 
            figure;
            clear Fx Fz
            
            %% ML force xcorr ML clavicle disp
                    
            subplot(numrows,numcols,1)
%             % Remove means
            Fx = TrialData(i).Results.Forces(1,3:end)-mean(TrialData(i).Results.Forces(1,3:end));
            swayX = TrialData(i).Results.CLAV(1,3:end)-nanmean(TrialData(i).Results.CLAV(1,3:end));
%             % Compute xcorr, get max value.
%             [r, lags] = xcorr(Fx,swayX,fs,'normalized');
%             ind = find(abs(r) == max(abs(r)));
%             if length(ind) > 1 
%                 ind = ind
%             end
%             lagMax = lags(ind);
            % Do polyfit between force and sway at lag that maximizes xcorr
            % to get y axis scaling factor for plot
            lagMax = TrialData(i).Results.xcorrX;
            maxRx = TrialData(i).Results.xcorrX;
            if lagMax < 0
                timeShift = TrialData(i).Results.time(3:end+lagMax);
                swayXShift = swayX(-lagMax+1:end);
                FxCorr = Fx(1:end+lagMax);
            elseif lagMax > 0
                timeShift = TrialData(i).Results.time(3+lagMax:end);
                swayXShift = swayX(1:end-lagMax);
                FxCorr = Fx(1+lagMax:end);
            else
                timeShift = TrialData(i).Results.time(3:end);
                swayXShift = swayX;
                FxCorr = Fx;
            end
            
            p = polyfit(swayXShift,FxCorr,1);
            
            % Plot force and sway
            plot(TrialData(i).Results.time(3:end),Fx,'k-',timeShift,polyval(p,swayXShift),'b-'),hold on; % Plot with mean subtracted

            % Labels and formatting
            ylabel('F ML (N)'); box off; set(gca,'tickdir','out');
            titlename = sprintf('%s %s: r=%.2f lag = %.2f',TrialData(i).Info.Condition,TrialData(i).Info.Trial,maxRx,lagMax/fs);
            legend('F','ML sway','orientation','horizontal','location','northoutside');
            title(titlename);
            
%             %% ML force, vertical clavicle displacement
%             
%             subplot(numrows,numcols,3)
%             % Remove mean
%             swayZ = TrialData(i).Results.CLAV(3,3:end)-nanmean(TrialData(i).Results.CLAV(3,3:end));
%             % Compute xcorr, get max value.
%             [r, lags] = xcorr(Fx,swayZ,fs,'normalized');
%             ind = find(abs(r) == max(abs(r)));
%             if length(ind) > 1 
%                 ind = ind
%             end
%             lagMax = lags(ind);
%             % Do polyfit between force and sway at lag that maximizes xcorr
%             % to get y axis scaling factor for plot
%             if lagMax < 0
%                 timeShift = TrialData(i).Results.time(3:end+lagMax);
%                 swayZShift = swayZ(-lagMax+1:end);
%                 FxCorr = Fx(1:end+lagMax);
%             elseif lagMax > 0
%                 timeShift = TrialData(i).Results.time(3+lagMax:end);
%                 swayZShift = swayZ(1:end-lagMax);
%                 FxCorr = Fx(1+lagMax:end);
%             else
%                 timeShift = TrialData(i).Results.time(3:end);
%                 swayZShift = swayZ;
%                 FxCorr = Fx;
%             end
%             
%             p = polyfit(swayZShift,FxCorr,1);
%             
%             % Plot force and sway
%             plot(TrialData(i).Results.time(3:end),Fx,'k-',timeShift,polyval(p,swayZShift),'b-'),hold on; % Plot with mean subtracted
% 
%             % Labels and formatting
%             ylabel('F ML (N)'); box off; set(gca,'tickdir','out');
%             titlename = sprintf('r=%.2f lag = %.2f',r(ind),lagMax/fs);
%             legend('F','Vertical sway','orientation','horizontal','location','northoutside');
%             title(titlename);
%                         
%             %% Vertical force xcorr ML clavicle disp
%                     
%             subplot(numrows,numcols,2)
%             % Remove means
%             Fz = TrialData(i).Results.Forces(3,3:end)-mean(TrialData(i).Results.Forces(3,3:end));
%             swayX = TrialData(i).Results.CLAV(1,3:end)-nanmean(TrialData(i).Results.CLAV(1,3:end));
%             % Compute xcorr, get max value.
%             [r, lags] = xcorr(Fz,swayX,fs,'normalized');
%             ind = find(abs(r) == max(abs(r)));
%             if length(ind) > 1 
%                 ind = ind
%             end
%             lagMax = lags(ind);
%             % Do polyfit between force and sway at lag that maximizes xcorr
%             % to get y axis scaling factor for plot
%             if lagMax < 0
%                 timeShift = TrialData(i).Results.time(3:end+lagMax);
%                 swayXShift = swayX(-lagMax+1:end);
%                 FzCorr = Fz(1:end+lagMax);
%             elseif lagMax > 0
%                 timeShift = TrialData(i).Results.time(3+lagMax:end);
%                 swayXShift = swayX(1:end-lagMax);
%                 FzCorr = Fz(1+lagMax:end);
%             else
%                 timeShift = TrialData(i).Results.time(3:end);
%                 swayXShift = swayX;
%                 FzCorr = Fz;
%             end
% 
%             p = polyfit(swayXShift,FzCorr,1);
%             
%             % Plot force and sway
%             plot(TrialData(i).Results.time(3:end),Fz,'k-',timeShift,polyval(p,swayXShift),'b-'),hold on; % Plot with mean subtracted
% 
%             % Labels and formatting
%             ylabel('F Vert (N)'); box off; set(gca,'tickdir','out');
%             titlename = sprintf('r=%.2f lag = %.2f',r(ind),lagMax/fs);
%             legend('F','ML sway','orientation','horizontal','location','northoutside');
%             title(titlename);
%             
            %% Vertical force xcorr Vert clavicle disp
                    
            subplot(numrows,numcols,4)
%             % Remove means
            Fz = TrialData(i).Results.Forces(3,3:end)-mean(TrialData(i).Results.Forces(3,3:end));
            swayZ = TrialData(i).Results.CLAV(3,3:end)-nanmean(TrialData(i).Results.CLAV(3,3:end));
%             % Compute xcorr, get max value.
%             [r, lags] = xcorr(Fz,swayZ,fs,'normalized');
%             ind = find(abs(r) == max(abs(r)));
%             if length(ind) > 1 
%                 ind = ind
%             end
%             lagMax = lags(ind);
            lagMax = TrialData(i).Results.xcorrZ;
            maxRz = TrialData(i).Results.xcorrZ;
            % Do polyfit between force and sway at lag that maximizes xcorr
            % to get y axis scaling factor for plot
            if lagMax < 0
                timeShift = TrialData(i).Results.time(3:end+lagMax);
                swayZShift = swayZ(-lagMax+1:end);
                FzCorr = Fz(1:end+lagMax);
            elseif lagMax > 0
                timeShift = TrialData(i).Results.time(3+lagMax:end);
                swayZShift = swayZ(1:end-lagMax);
                FzCorr = Fz(1+lagMax:end);
            else
                timeShift = TrialData(i).Results.time(3:end);
                swayZShift = swayZ;
                FzCorr = Fz;
            end
            
            p = polyfit(swayZShift,FzCorr,1);
            
            % Plot force and sway
            plot(TrialData(i).Results.time(3:end),Fz,'k-',timeShift,polyval(p,swayZShift),'b-'),hold on; % Plot with mean subtracted

            % Labels and formatting
            ylabel('F Vert (N)'); box off; set(gca,'tickdir','out');
            titlename = sprintf('r=%.2f lag = %.2f',maxRz,lagMax/fs);
            legend('F','Vert sway','orientation','horizontal','location','northoutside');
            title(titlename);
            
            %% Save one pdf per subject
            set(gcf,'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
            pdfname = sprintf('HHI%i_force_clav',subj);
            export_fig(gcf,pdfname,'-pdf','-append')
            close all;
        end
    end
end
