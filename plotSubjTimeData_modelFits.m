% Plot time series data for an individual subject. Plot only significant
% predictors from regression. Tension is positive for sensor in all
% directions. However, need to negate force in all directions when plotting
% force on POB.

% Use this code to plot example subject data for model fits
% Fit F to torso state: HHI08 t15 and HHI
% Fit F to IP state: HHI11 t23 and HHI10 t2
clear; clc; %close all;

subj_array_force = [3:5 8:13];
s2 = [5 1 11 8.5]; % paper size
plotind = 1;
numcols = 2; numrows = 4;
colors(1,:) = [0.00,0.45,0.74]; % nice blue
colors(2,:) = [0.85,0.33,0.10]; % nice red
colors(3,:) = [0.47,0.67,0.19]; % nice green

plotIP = 1; % regress to IP state

%%
for subj = 11%subj_array_force
    filename = sprintf('HHI2017_%i.mat',subj);
    load(filename);

    % Cycle through all trials and plot time series data for actual force,
    % predicted force, marker data, etc.

    %% ML/x dir L col

    for i = 23%1:length(TrialData)
       
        if strcmp(TrialData(i).Info.Condition,'Assist Beam') %strcmp(TrialData(i).Info.Condition,'Solo Ground') %|| strcmp(TrialData(i).Info.Condition,'Assist Ground') 
%             figure;
            
            %% ML dir
                    
%             plotind = 1;            
            
            if plotIP == 1
                % Plot actual force and components
                subplot(numrows,numcols,plotind)
                clear temp m;
                temp.c = TrialData(i).Results.cx_IP;
                temp.c(isnan(temp.c)) = 0;
                m = [TrialData(i).Results.IntPtAcc(:,1) TrialData(i).Results.IntPtVel(2:end,1) TrialData(i).Results.IntPt(3:end,1)-nanmean(TrialData(i).Results.IntPt(3:end,1))];
                tAnalysis = TrialData(i).Results.time(3:end); % only can do regression for 3:end since need acc
                plot(tAnalysis,TrialData(i).Results.Forces(3:end,1)-temp.c(4),'color',0.75.*[1 1 1],'linewidth',2);hold on
                plot(tAnalysis,m*temp.c(1:3),'k');
                plot(tAnalysis,TrialData(i).Results.cx_IP(1).*TrialData(i).Results.IntPtAcc(:,1),'color',colors(1,:));
                plot(tAnalysis,TrialData(i).Results.cx_IP(2).*TrialData(i).Results.IntPtVel(2:end,1),'color',colors(2,:));
                plot(tAnalysis,TrialData(i).Results.cx_IP(3).*(TrialData(i).Results.IntPt(3:end,1)-nanmean(TrialData(i).Results.IntPt(3:end,1))),'color',colors(3,:));
                % Labels and formatting
                if plotind == 1
                    legend('F','Fmodel','Facc','Fvel','Fdisp','orientation','horizontal');
                end
                titlename = sprintf('HHI%i %s R^2 = %.2f, m = %.2f, b = %.2f, k = %.2f',subj,TrialData(i).Info.Trial,TrialData(i).Results.rsqx_torso,TrialData(i).Results.cx_torso(1),TrialData(i).Results.cx_torso(2),TrialData(i).Results.cx_torso(3));
                title(titlename);
                ylabel('Force (N)'); box off; set(gca,'tickdir','out');
                xlabel('Time (s)');

                % Plot IP acc
                plotind = plotind + numcols; 
                subplot(numrows,numcols,plotind),hold on;
                plot(tAnalysis,TrialData(i).Results.IntPtAcc(:,1),'color',colors(1,:))
                hline(0,'k--');
                % Labels and formatting
                ylabel('xIddot (m/s^2)'); box off; set(gca,'tickdir','out');

                % Plot IP vel
                plotind = plotind + numcols; 
                subplot(numrows,numcols,plotind),hold on;
                plot(tAnalysis,TrialData(i).Results.IntPtVel(2:end,1),'color',colors(2,:))
                hline(0,'k--');
                % Labels and formatting
                ylabel('xIdot (m/s)'); box off; set(gca,'tickdir','out');

                % Plot IP disp
                plotind = plotind + numcols; 
                subplot(numrows,numcols,plotind),hold on;
                plot(tAnalysis,TrialData(i).Results.IntPt(3:end,1)-nanmean(TrialData(i).Results.IntPt(3:end,1)),'color',colors(3,:))
                hline(0,'k--');
                % Labels and formatting
                ylabel('xI (m)'); box off; set(gca,'tickdir','out');
            end

        end
    end
end