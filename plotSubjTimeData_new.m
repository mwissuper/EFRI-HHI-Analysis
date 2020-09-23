 % Plot time series for all trials per subj for given variable to check metrics
% were calculated correctly. Plot all trials, including early trials

clear all; clc; close all;

subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)
plotDistSpeed = 0;
plotSway = 0;
plotIPstate = 0; % plot IP pos and vel to check if marker gapfill error or other marker error. Used in power calc's
plotFvIP = 0; % Plot IV vel and force to check signed force metrics and power metrics
plotFmag = 0;
plotF = 0; % signed force
plotV = 0; % signed vel IP
plotPmag = 0;
plotP = 0; % signed power
plotStrategy = 0; % Opp/amp dev/ret
plotTorsoFit = 1; % plot without the constant term
plotIPFit = 0; % plot without the constant term

colors(1,:) = [0.00,0.45,0.74]; % nice blue
colors(2,:) = [0.85,0.33,0.10]; % nice red
colors(3,:) = [0.47,0.67,0.19]; % nice green

if plotSway == 1 % 3 cond's
    numrows = 5; numcols = 6;
elseif plotIPstate == 1 || plotFvIP == 1 || plotFmag == 1 || plotPmag == 1 || plotStrategy == 1 || plotV == 1 || plotP == 1 || plotF == 1% 2 cond's
    numrows = 4; numcols = 5;
elseif plotTorsoFit == 1 || plotIPFit == 1 % Assist Beam cond only
    numrows = 4; numcols = 3;
end

%% Cycle through all subj's and plot all trials
% Plot percent of trial with combo's of signs in F and P vs. trial

subjind = 0;
temp = [];
for subj = 3% subj_array
    clear Fx Fy Fz
    subjind = subjind + 1;
    filename = sprintf('HHI2017_%i.mat',subj);
    load(filename);
    plotind = 0;
    sample_rate = TrialData(1).Markers.samplerate;
    % get beamMidline from mean of all indiv beam trials' midlines
    temp = [];
    for i = 1:length(TrialData)
        if any(strcmp(TrialData(i).Info.Condition,{'Solo Beam','Assist Beam'}))
            a = TrialData(i).Results.beamMidline/1000;
            temp = [temp a];
        end
    end
    beamMidline = nanmean(temp);
    for i = 1:length(TrialData)
        if ~strcmp(TrialData(i).Info.Condition,'Assist Solo') & ~strcmp(TrialData(i).Info.Condition,'FH') & ~strcmp(TrialData(i).Info.Condition,'AP') & ~strcmp(TrialData(i).Info.Condition,'POB') & ~strcmp(TrialData(i).Info.Condition,'Force Handle')  % disregard these trials
            if strcmp(TrialData(i).Info.Condition,'Solo Beam') || strcmp(TrialData(i).Info.Condition,'Solo Ground')
                clav = (TrialData(i).Markers.CLAV-TrialData(i).Markers.CLAV(TrialData(i).Results.startIdx,:))/1000;
            else
                clav = (TrialData(i).Markers.POB.CLAV-TrialData(i).Markers.POB.CLAV(TrialData(i).Results.startIdx,:))/1000;
                IP = (TrialData(i).Markers.AP.RFIN-TrialData(i).Markers.AP.RFIN(TrialData(i).Results.startIdx,:))/1000;
            end
            t = 0:(length(clav)-1);
            t = t./sample_rate;
            tAnalysis = TrialData(i).Results.time + t(TrialData(i).Results.startIdx); % Shift results time so match marker data
            if plotDistSpeed == 1
                if any(strcmp(TrialData(i).Info.Condition,{'Solo Beam','Assist Beam'}))
                    plotind = plotind + 1; hold on;
                    tn = str2num(TrialData(i).Info.Trial(end-1:end));
                    yyaxis left
                    if strcmp(TrialData(i).Info.Condition,'Solo Beam')
                        plot(tn,TrialData(i).Results.totalDistance,'x')
                    elseif strcmp(TrialData(i).Info.Condition,'Assist Beam')
                        plot(tn,TrialData(i).Results.totalDistance,'o')
                    end
                    if plotind == 1
                       ylabel('Distance completed (m)'),xlabel('Trial');
                       set(gcf,'outerposition',[672   662   576   389]);
                       titlename = sprintf('HHI%i',subj); title(titlename); 
                    end
                    yyaxis right
                    if strcmp(TrialData(i).Info.Condition,'Solo Beam')
                        plot(tn,TrialData(i).Results.avgSpeed,'x')
                    elseif strcmp(TrialData(i).Info.Condition,'Assist Beam')
                        plot(tn,TrialData(i).Results.avgSpeed,'o')
                    end
                    if plotind == 1
                       ylabel('Avg speed (m/s)'); 
                    end 
                end
            elseif plotSway == 1
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                plot(t,clav(:,1)),hold on;
                plot(tAnalysis,TrialData(i).Results.beamerSway); % could be clav or c7
                plot(tAnalysis,(mean(TrialData(i).Results.beamerSway)+std(TrialData(i).Results.beamerSway))*ones(size(tAnalysis))) % Plot line for SD
                plot(tAnalysis,(mean(TrialData(i).Results.beamerSway)-std(TrialData(i).Results.beamerSway))*ones(size(tAnalysis))) 
                xlabel('Time (s)'),ylabel('Clav sway (m)');
                if plotind == 1
                    legend('clav','torso');
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
            elseif plotIPstate == 1 && ~strcmp(TrialData(i).Info.Condition,'Solo Beam')
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                yyaxis left
                plot(t,IP(:,1)),hold on;
                plot(tAnalysis,TrialData(i).Results.IntPt(:,1)-TrialData(i).Results.IntPt(1,1),'g-');
                xlabel('Time (s)'),ylabel('Lateral IP pos (m)');
                yyaxis right
                plot(tAnalysis(2:end),TrialData(i).Results.IntPtVel(:,1));
                ylabel('Lateral IP vel (m/s)');
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                    legend('raw marker','used for calc','filt','orientation','horizontal');
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
            elseif plotFvIP == 1 && ~strcmp(TrialData(i).Info.Condition,'Solo Beam')
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                yyaxis left
                plot(tAnalysis(2:end),TrialData(i).Results.IntPtVel(:,1));,hold on
                hline(0,'k-')
                if find(mod(plotind,numcols) == 1)
                    ylabel('Lateral IP vel (m/s)');
                end
                if subj == 9
                    ylim([-0.4 0.4]);
                end
                yyaxis right
                plot(tAnalysis,TrialData(i).Results.Forces(:,1));
                if subj == 9
                    ylim([-20 20]);
                end
                if find(mod(plotind,numcols) == 0)
                    ylabel('Lateral IP F (N)');
                end
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
            elseif plotFmag == 1 && ~strcmp(TrialData(i).Info.Condition,'Solo Beam')
                if subj == 12
                    numrows = 5; numcols = 5;
                else
                    numrows = 4; numcols = 5;
                end
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                plot(tAnalysis,abs(TrialData(i).Results.Forces(:,1)));
                hline(nanmean(abs(TrialData(i).Results.Forces(:,1))),'k--');
                ylim([0 20]);
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
                if mod(plotind,numcols) == 1
                    ylabel('Force (N)')
                end
                xlabel('Time (s)');
            elseif plotF == 1 && ~strcmp(TrialData(i).Info.Condition,'Solo Beam') && ~strcmp(TrialData(i).Info.Condition,'Solo Ground')
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                plot(tAnalysis,TrialData(i).Results.Forces(:,1));
                hline(nanmean(TrialData(i).Results.Forces(:,1)),'k--');
                xlabel('Time (s)'),ylabel('Force (N)');
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
            elseif plotV == 1 && ~strcmp(TrialData(i).Info.Condition,'Solo Beam') && ~strcmp(TrialData(i).Info.Condition,'Solo Ground')
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                plot(tAnalysis(2:end),TrialData(i).Results.IntPtVel(:,1));
                xlabel('Time (s)'),ylabel('Lateral IP vel (m/s)');
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
            elseif plotP == 1 && ~strcmp(TrialData(i).Info.Condition,'Solo Beam') && ~strcmp(TrialData(i).Info.Condition,'Solo Ground')
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                plot(tAnalysis(2:end),TrialData(i).Results.IntPower(:,1));
                hline(nanmean(TrialData(i).Results.IntPower(:,1)),'k--');
                xlabel('Time (s)'),ylabel('Power (W)');
    %                 if subj == 10
%                     ylim([0 4]);
    %                 end
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
             elseif plotStrategy == 1 && (strcmp(TrialData(i).Info.Condition,'Assist Ground') || strcmp(TrialData(i).Info.Condition,'Assist Beam'))
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                yyaxis left
                if strcmp(TrialData(i).Info.Condition,'Assist Ground')
                    midline = nanmean(TrialData(i).Results.torso(:,1));
                elseif strcmp(TrialData(i).Info.Condition,'Assist Beam')
                    midline = beamMidline;
                end
                plot(tAnalysis(2:end),TrialData(i).Results.torso(2:end,1)-midline);%.*TrialData(i).Results.vTorso(:,1));
%                 hline(temp.midline,'b--');
                ylim([-.1 .1]); % torso disp
                xlabel('Time (s)'),ylabel('Torso disp (m)') %ylabel('x_b*dx_b (m^2/s)'),ylim([-.2 .2])
                yyaxis right
                if strcmp(TrialData(i).Info.Condition,'Assist Ground') || strcmp(TrialData(i).Info.Condition,'Assist Beam') 
                    if ~isnan(TrialData(i).Results.Forces) % good trials
                        plot(tAnalysis(2:end),TrialData(i).Results.Forces(2:end,1));%.*TrialData(i).Results.vTorso(:,1));
    %                     ylabel('Power torso (W)')
                        ylabel('IP force (N)')
                    end
                end
                if subj == 5
                    ylim([-10 10]); % Force
                    hline(0,'r--');
                end
%                 axis tight
%                 ylim([-1.5 1.5]); % power
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
            elseif plotTorsoFit == 1 & strcmp(TrialData(i).Info.Condition,'Assist Beam') & (isnan(TrialData(i).Results.IntPower) == 0)
                plotind = plotind + 1;
                tAnalysis(1:2) = []; % only can do regression for 3:end since need acc
                clear temp;
                temp.c = TrialData(i).Results.cx_torso;
                temp.c(isnan(temp.c)) = 0;
                m = [TrialData(i).Results.aTorso(:,1) TrialData(i).Results.vTorso(2:end,1) TrialData(i).Results.torso(3:end,1)-TrialData(i).Results.beamMidline];
                subplot(numrows,numcols,plotind)
                plot(tAnalysis,TrialData(i).Results.Forces(3:end,1)-temp.c(4),'color',0.75.*[1 1 1],'linewidth',2);hold on
                plot(tAnalysis,m*temp.c(1:3),'k');
                plot(tAnalysis,TrialData(i).Results.cx_torso(1).*TrialData(i).Results.aTorso(:,1),'color',colors(1,:));
                plot(tAnalysis,TrialData(i).Results.cx_torso(2).*TrialData(i).Results.vTorso(2:end,1),'color',colors(2,:));
                plot(tAnalysis,TrialData(i).Results.cx_torso(3).*(TrialData(i).Results.torso(3:end,1)-TrialData(i).Results.beamMidline),'color',colors(3,:));
                if plotind == 1
                    titlename = sprintf('HHI%i %s R^2 = %.2f',subj,TrialData(i).Info.Trial,TrialData(i).Results.rsqx_torso);
                    legend('F','Fmodel','Facc','Fvel','Fdisp','orientation','horizontal');
                else
                    titlename = sprintf('%s R^2 = %.2f, m = %.2f, b = %.2f, k = %.2f',TrialData(i).Info.Trial,TrialData(i).Results.rsqx_torso,TrialData(i).Results.cx_torso(1),TrialData(i).Results.cx_torso(2),TrialData(i).Results.cx_torso(3));
                end
                title(titlename);xlabel('Time (s)'),ylabel('Lateral Force (N)');
            elseif plotIPFit == 1 & strcmp(TrialData(i).Info.Condition,'Assist Beam') & (isnan(TrialData(i).Results.IntPower) == 0)
                plotind = plotind + 1;
                tAnalysis(1:2) = []; % only can do regression for 3:end since need acc
                clear temp;
                temp.c = TrialData(i).Results.cx_torso;
                temp.c(isnan(temp.c)) = 0;
                m = [TrialData(i).Results.IntPtAcc(:,1) TrialData(i).Results.IntPtVel(2:end,1) TrialData(i).Results.IntPt(3:end,1)-nanmean(TrialData(i).Results.IntPt(3:end,1))];
                subplot(numrows,numcols,plotind)
                plot(tAnalysis,TrialData(i).Results.Forces(3:end,1)-temp.c(4),'color',0.75.*[1 1 1],'linewidth',2);hold on
                plot(tAnalysis,m*temp.c(1:3),'k');
                plot(tAnalysis,TrialData(i).Results.cx_IP(1).*TrialData(i).Results.IntPtAcc(:,1),'color',colors(1,:));
                plot(tAnalysis,TrialData(i).Results.cx_IP(2).*TrialData(i).Results.IntPtVel(2:end,1),'color',colors(2,:));
                plot(tAnalysis,TrialData(i).Results.cx_IP(3).*(TrialData(i).Results.IntPt(3:end,1)-nanmean(TrialData(i).Results.IntPt(3:end,1))),'color',colors(3,:));
                if plotind == 1
                    titlename = sprintf('HHI%i %s R^2 = %.2f, m = %.2f, b = %.2f, k = %.2f',subj,TrialData(i).Info.Trial,TrialData(i).Results.rsqx_IP,TrialData(i).Results.cx_IP(1),TrialData(i).Results.cx_IP(2),TrialData(i).Results.cx_IP(3));
                    legend('F','Fmodel','Facc','Fvel','Fdisp','orientation','horizontal');
                else
                    titlename = sprintf('%s R^2 = %.2f, m = %.2f, b = %.2f, k = %.2f',TrialData(i).Info.Trial,TrialData(i).Results.rsqx_IP,TrialData(i).Results.cx_IP(1),TrialData(i).Results.cx_IP(2),TrialData(i).Results.cx_IP(3));
                end
                title(titlename);xlabel('Time (s)'),ylabel('Lateral Force (N)');
            end
        end
    end
    
    a = findobj(gcf,'type','axes');
    set(a,'tickdir','out','box','off')

    if plotDistSpeed == 1
        figname = sprintf('HHI%i_dist_speed',subj);
    elseif plotSway == 1
        figname = sprintf('HHI%i_sway',subj);
    elseif plotIPstate == 1
        figname = sprintf('HHI%i_IPstate',subj);
    elseif plotFmag == 1
        figname = sprintf('HHI%i_FmagVar',subj);
    elseif plotFvIP == 1
        figname = sprintf('HHI%i_FvIP',subj);
    elseif plotTorsoFit == 1
        figname = sprintf('HHI%i_forceFitTorso',subj);
    end
    if plotDistSpeed ~= 1
        set(gcf,'outerposition',[  -7         180        1072         909]);
    end
    saveas(gcf,figname,'fig');
    close all;
end
