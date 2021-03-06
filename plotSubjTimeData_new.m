 % Plot time series for all trials per subj for given variable to check metrics
% were calculated correctly. Plot all trials, including early trials

clear; clc; close all;

subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)
plotDistSpeed = 0;
plotSway = 0;
plotIPstate = 0; % plot IP pos and vel to check if marker gapfill error or other marker error. Used in power calc's
plotFvIP = 0; % Plot IV vel and force to check signed force metrics and power metrics
plotFmag = 0;
plotF = 0; % signed force
plotT = 0; % signed Ty and torso state vectors to check signs make sense
plotTFxFz = 1; % Plot signed Ty due to Fx and Fz components separately
plotV = 0; % signed vel IP
plotPmag = 0;
plotP = 0; % signed power
plotPang = 0; % Angular power (signed)
plotPangRMS = 0; % Angular power (RMS)
plotStrategy = 0; % Opp/amp dev/ret
plotTorsoFitX = 0; % plot without the constant term
plotTorsoFitZ = 0; % plot without the constant term
plotTorsoFitAng = 0; % plot without the constant term
titleVAF = 0; % Change title to display VAF of each term of model (1) instead of value of coeff (0) or Rsq of components (2)
plotIPFit = 0; % plot without the constant term

colors(1,:) = [0.00,0.45,0.74]; % nice blue
colors(2,:) = [0.85,0.33,0.10]; % nice red
colors(3,:) = [0.47,0.67,0.19]; % nice green

if plotSway == 1 % 3 cond's
    numrows = 5; numcols = 6;
elseif plotIPstate == 1 || plotFvIP == 1 || plotFmag == 1 || plotPmag == 1 || plotStrategy == 1 || plotV == 1 || plotP == 1 || plotF == 1 || plotT == 1 || plotPang == 1 || plotPangRMS == 1 || plotTFxFz == 1 % 2 cond's
    numrows = 4; numcols = 5;
elseif plotTorsoFitX == 1 || plotTorsoFitZ == 1 || plotIPFit == 1  || plotTorsoFitAng == 1 %|| plotTorsoFitAng == 1 % Assist Beam cond only
    numrows = 4; numcols = 3;
end

%% Cycle through all subj's and plot all trials
% Plot percent of trial with combo's of signs in F and P vs. trial

subjind = 0;
temp = [];
for subj = subj_array_force
    if subj == 12
        numrows = numrows + 1; 
    end
    clear Fx Fy Fz 
    VAF = []; Rsq = [];
    subjind = subjind + 1;
    filename = sprintf('HHI2017_%i.mat',subj);
    load(filename);
    plotind = 0;
    sample_rate = TrialData(1).Markers.samplerate;

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
                plot(tAnalysis,TrialData(i).Results.Force(:,1));
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
                plot(tAnalysis,abs(TrialData(i).Results.Force(:,1)));
                hline(nanmean(abs(TrialData(i).Results.Force(:,1))),'k--');
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
                plot(tAnalysis,TrialData(i).Results.Force(:,1));
                hline(nanmean(TrialData(i).Results.Force(:,1)),'k--');
                xlabel('Time (s)'),ylabel('Force (N)');
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
            elseif plotT == 1 && ~strcmp(TrialData(i).Info.Condition,'Solo Beam') && ~strcmp(TrialData(i).Info.Condition,'Solo Ground') % Plot T total, Tg, Text
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                plot(tAnalysis,TrialData(i).Results.TyGd,tAnalysis,TrialData(i).Results.Tg,tAnalysis,TrialData(i).Results.TyExt),hold on;
                xlabel('Time (s)'),ylabel('Ty gd (Nm)');
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                    legend('Total','Gravity','Ext');
                else
                    titlename = sprintf('%s %s',TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                end
                title(titlename);
            elseif plotTFxFz == 1 & ~strcmp(TrialData(i).Info.Condition,'Solo Beam') & ~strcmp(TrialData(i).Info.Condition,'Solo Ground') & ~isnan(TrialData(i).Results.Force)
                clear r TFx TFz
                
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                
                r = TrialData(i).Results.IP - [TrialData(i).Results.midline 0 0];
                TFx = r(:,3).*TrialData(i).Results.Force(:,1);
                TFz = -r(:,1).*TrialData(i).Results.Force(:,3);

                plot(tAnalysis,TFx,tAnalysis,TFz,tAnalysis,TrialData(i).Results.TyGd),hold on;
                ylim([-10 10]); hline(0,'k-');
                xlabel('Time (s)'),ylabel('Ty about gd/beam (Nm)');
                
                if plotind == 1
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(i).Info.Trial,TrialData(i).Info.Condition);
                    legend('Fx comp','Fz comp','Total','orientation','horizontal');
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
            elseif (plotPang == 1 || plotPangRMS == 1) && ~strcmp(TrialData(i).Info.Condition,'Solo Beam') && ~strcmp(TrialData(i).Info.Condition,'Solo Ground')
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                if plotPang == 1
                    plot(tAnalysis,TrialData(i).Results.angP);
                    hline(nanmean(TrialData(i).Results.angP),'k--');
                    ylim([-0.1 0.1]); % Make ylim's small to see diff in means
                else % RMS
                    plot(tAnalysis,sqrt(TrialData(i).Results.angP.^2));
                    hline(nanmean(sqrt(TrialData(i).Results.angP.^2)),'k--');
                end
                xlabel('Time (s)'),ylabel('Power (N*m/s)');
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
                    if ~isnan(TrialData(i).Results.Force) % good trials
                        plot(tAnalysis(2:end),TrialData(i).Results.Force(2:end,1));%.*TrialData(i).Results.vTorso(:,1));
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
            elseif plotTorsoFitX == 1 & strcmp(TrialData(i).Info.Condition,'Assist Beam') & (isnan(TrialData(i).Results.IntPower) == 0)
                plotind = plotind + 1;
                tAnalysis(1:2) = []; % only can do regression for 3:end since need acc
                Fa = TrialData(i).Results.cx_torso(1).*TrialData(i).Results.aTorso(:,1);
                Fv = TrialData(i).Results.cx_torso(2).*TrialData(i).Results.vTorso(2:end,1);
                Fp = TrialData(i).Results.cx_torso(3).*(TrialData(i).Results.torso(3:end,1)-TrialData(i).Results.beamMidline);
                clear temp;
                temp.c = TrialData(i).Results.cx_torso;
                temp.c(isnan(temp.c)) = 0;
                m = [TrialData(i).Results.aTorso(:,1) TrialData(i).Results.vTorso(2:end,1) TrialData(i).Results.torso(3:end,1)-TrialData(i).Results.beamMidline];
                subplot(numrows,numcols,plotind)
                plot(tAnalysis,TrialData(i).Results.Force(3:end,1)-temp.c(4),'color',0.75.*[1 1 1],'linewidth',2);hold on
                plot(tAnalysis,m*temp.c(1:3),'k');
                plot(tAnalysis,Fa,'color',colors(1,:));
                plot(tAnalysis,Fv,'color',colors(2,:));
                plot(tAnalysis,Fp,'color',colors(3,:));
                % Calculate VAF
                VAF_m = rsqr_uncentered(TrialData(i).Results.Force(3:end,1)-temp.c(4),Fa);
                VAF_b = rsqr_uncentered(TrialData(i).Results.Force(3:end,1)-temp.c(4),Fv);
                VAF_k = rsqr_uncentered(TrialData(i).Results.Force(3:end,1)-temp.c(4),Fp);
                % Calculate Rsq
                Rsq_m = rsqr(TrialData(i).Results.Force(3:end,1)-temp.c(4),Fa);
                Rsq_b = rsqr(TrialData(i).Results.Force(3:end,1)-temp.c(4),Fv);
                Rsq_k = rsqr(TrialData(i).Results.Force(3:end,1)-temp.c(4),Fp);
                % Save these values for each subj
                VAF = [VAF; VAF_m VAF_b VAF_k];
                Rsq = [Rsq; Rsq_m Rsq_b Rsq_k];
                
                if plotind == 1
                    if titleVAF == 1
                        titlename = sprintf('HHI%i %s, VAF: m = %.2f, b = %.2f, k = %.2f',subj,TrialData(i).Info.Trial,VAF_m,VAF_b,VAF_k);
                    elseif titleVAF == 2
                        titlename = sprintf('HHI%i %s, R^2: %m = %.2f, b = %.2f, k = %.2f',subj,TrialData(i).Info.Trial,Rsq_m,Rsq_b,Rsq_k);
                    else
                        titlename = sprintf('HHI%i %s R^2 = %.2f, m = %.2f, b = %.2f, k = %.2f',subj,TrialData(i).Info.Trial,TrialData(i).Results.rsqx_torso,TrialData(i).Results.cx_torso(1),TrialData(i).Results.cx_torso(2),TrialData(i).Results.cx_torso(3));
                    end
                    legend('F','Fmodel','Facc','Fvel','Fdisp','orientation','horizontal');
                else
                    if titleVAF == 1
                        titlename = sprintf('%s , VAF: m = %.2f, b = %.2f, k = %.2f',TrialData(i).Info.Trial,VAF_m,VAF_b,VAF_k);
                    elseif titleVAF == 2
                        titlename = sprintf('%s R^2: m = %.2f, b = %.2f, k = %.2f',TrialData(i).Info.Trial,Rsq_m,Rsq_b,Rsq_k);
                    else
                        titlename = sprintf('%s R^2 = %.2f, m = %.2f, b = %.2f, k = %.2f',TrialData(i).Info.Trial,TrialData(i).Results.rsqx_torso,TrialData(i).Results.cx_torso(1),TrialData(i).Results.cx_torso(2),TrialData(i).Results.cx_torso(3));
                    end
                end
                title(titlename);xlabel('Time (s)'),ylabel('Lateral Force (N)');
            elseif plotTorsoFitZ == 1 & strcmp(TrialData(i).Info.Condition,'Assist Beam') & (isnan(TrialData(i).Results.IntPower) == 0)
                plotind = plotind + 1;
                tAnalysis(1:2) = []; % only can do regression for 3:end since need acc
                Fa = TrialData(i).Results.cz_torso(1).*TrialData(i).Results.aTorso(:,3);
                Fv = TrialData(i).Results.cz_torso(2).*TrialData(i).Results.vTorso(2:end,3);
                Fp = TrialData(i).Results.cz_torso(3).*(TrialData(i).Results.torso(3:end,3)-TrialData(i).Results.eqbmHt);
                temp = [];
                temp.c = TrialData(i).Results.cz_torso;
                temp.c(isnan(temp.c)) = 0;
                m = [TrialData(i).Results.aTorso(:,3) TrialData(i).Results.vTorso(2:end,3) TrialData(i).Results.torso(3:end,3)-TrialData(i).Results.eqbmHt];
                subplot(numrows,numcols,plotind)
                plot(tAnalysis,TrialData(i).Results.Force(3:end,3)-temp.c(4),'color',0.75.*[1 1 1],'linewidth',2);hold on
                plot(tAnalysis,m*temp.c(1:3),'k');
                plot(tAnalysis,Fa,'color',colors(1,:));
                plot(tAnalysis,Fv,'color',colors(2,:));
                plot(tAnalysis,Fp,'color',colors(3,:));
                % Calculate VAF
                VAF_m = rsqr_uncentered(TrialData(i).Results.Force(3:end,3)-temp.c(4),Fa);
                VAF_b = rsqr_uncentered(TrialData(i).Results.Force(3:end,3)-temp.c(4),Fv);
                VAF_k = rsqr_uncentered(TrialData(i).Results.Force(3:end,3)-temp.c(4),Fp);
                % Calculate Rsq
                Rsq_m = rsqr(TrialData(i).Results.Force(3:end,3)-temp.c(4),Fa);
                Rsq_b = rsqr(TrialData(i).Results.Force(3:end,3)-temp.c(4),Fv);
                Rsq_k = rsqr(TrialData(i).Results.Force(3:end,3)-temp.c(4),Fp);
                % Save these values for each subj
                VAF = [VAF; VAF_m VAF_b VAF_k];
                Rsq = [Rsq; Rsq_m Rsq_b Rsq_k];
                
                if plotind == 1
                    if titleVAF == 1
                        titlename = sprintf('HHI%i %s, VAF: m = %.2f, b = %.2f, k = %.2f',subj,TrialData(i).Info.Trial,VAF_m,VAF_b,VAF_k);
                    elseif titleVAF == 2
                        titlename = sprintf('HHI%i %s, R^2: m = %.2f, b = %.2f, k = %.2f',subj,TrialData(i).Info.Trial,Rsq_m,Rsq_b,Rsq_k);
                    else
                        titlename = sprintf('HHI%i %s R^2 = %.2f, m = %.2f, b = %.2f, k = %.2f',subj,TrialData(i).Info.Trial,TrialData(i).Results.rsqz_torso,TrialData(i).Results.cz_torso(1),TrialData(i).Results.cz_torso(2),TrialData(i).Results.cz_torso(3));
                    end
                    legend('F','Fmodel','Facc','Fvel','Fdisp','orientation','horizontal');
                else
                    if titleVAF == 1
                        titlename = sprintf('%s , VAF: m = %.2f, b = %.2f, k = %.2f',TrialData(i).Info.Trial,VAF_m,VAF_b,VAF_k);
                    elseif titleVAF == 2
                        titlename = sprintf('%s R^2: m = %.2f, b = %.2f, k = %.2f',TrialData(i).Info.Trial,Rsq_m,Rsq_b,Rsq_k);
                    else
                        titlename = sprintf('%s R^2 = %.2f, m = %.2f, b = %.2f, k = %.2f',TrialData(i).Info.Trial,TrialData(i).Results.rsqz_torso,TrialData(i).Results.cz_torso(1),TrialData(i).Results.cz_torso(2),TrialData(i).Results.cz_torso(3));
                    end
                end
                title(titlename);xlabel('Time (s)'),ylabel('Vertical Force (N)');
            elseif plotTorsoFitAng == 1 & any(strcmpi(TrialData(i).Info.Condition,{'Assist Beam','Assist beam'})) & (isnan(TrialData(i).Results.angP) == 0)
                clear m temp
                ref = pi/2; 
%                 temp.theta = TrialData(i).Results.angTorso;
%                 indRev = find(temp.theta > ref);
%                 temp.alpha = TrialData(i).Results.alphaTorso;
%                 temp.alpha(indRev) = -TrialData(i).Results.alphaTorso(indRev);
%                 temp.w = TrialData(i).Results.wTorso;
%                 temp.w(indRev) = -TrialData(i).Results.wTorso(indRev);
                % No lag
                m = [TrialData(i).Results.alphaTorso TrialData(i).Results.wTorso TrialData(i).Results.angTorso-ref]; % Want displacement, not abs pos for fitting. For now, look at changes in pos of interaction point
%                 m = [temp.alpha temp.w temp.theta-ref]; % Want displacement, not abs pos for fitting. For now, look at changes in pos of interaction point
                temp.c = TrialData(i).Results.c_ang_torso;
                temp.c(isnan(temp.c)) = 0;
                
                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
%                 plot(tAnalysis,TrialData(i).Results.angTorso) % Angle should hover around pi/2
%                 hline(ref,'k-');
                plot(tAnalysis,TrialData(i).Results.TyGd-temp.c(4),'color',0.75.*[1 1 1],'linewidth',2);hold on
                plot(tAnalysis,m*temp.c(1:3),'k');
                plot(tAnalysis,TrialData(i).Results.c_ang_torso(1).*TrialData(i).Results.alphaTorso,'color',colors(1,:)); % Subtract out offset to better see quality of fits
                plot(tAnalysis,TrialData(i).Results.c_ang_torso(2).*TrialData(i).Results.wTorso,'color',colors(2,:));
                plot(tAnalysis,TrialData(i).Results.c_ang_torso(3).*(TrialData(i).Results.angTorso-ref),'color',colors(3,:));
%                 plot(tAnalysis,TrialData(i).Results.c_ang_torso(1).*temp.alpha,'color',colors(1,:)); % Subtract out offset to better see quality of fits
%                 plot(tAnalysis,TrialData(i).Results.c_ang_torso(2).*temp.w,'color',colors(2,:));
%                 plot(tAnalysis,TrialData(i).Results.c_ang_torso(3).*(temp.theta-ref),'color',colors(3,:));
                axis tight
                if plotind == 1
                    titlename = sprintf('HHI%i %s, R^2 = %.2f, I = %.2f, B = %.1f, K = %.0f',subj,TrialData(i).Info.Trial,TrialData(i).Results.rsq_ang_torso,TrialData(i).Results.c_ang_torso(1),TrialData(i).Results.c_ang_torso(2),TrialData(i).Results.c_ang_torso(3));
%                     legend('T','Tmodel','Talpha','Tw','Tang','orientation','horizontal');
                    legend('T ext','T model','T inertia','T damping','T stiffness','orientation','horizontal');
                else
                    titlename = sprintf('%s R^2 = %.2f, I = %.2f, B = %.1f, K = %.0f',TrialData(i).Info.Trial,TrialData(i).Results.rsq_ang_torso,TrialData(i).Results.c_ang_torso(1),TrialData(i).Results.c_ang_torso(2),TrialData(i).Results.c_ang_torso(3));
                end
                title(titlename);xlabel('Time (s)'),ylabel('Torque (Nm)');
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
    elseif plotTorsoFitX == 1
        figname = sprintf('HHI%i_forceFitTorso',subj);
        fname = sprintf('HHI%i_VAF_Rsq',subj);
    elseif plotTorsoFitAng == 1
        figname = sprintf('HHI%i_TorqueFitTorsoAng',subj);
    elseif plotTFxFz == 1
        figname = sprintf('HHI%i_TorqueFxFz',subj);
    end
    if plotDistSpeed ~= 1
        set(gcf,'outerposition',[  -7         172        1178         909]);
        set(gcf,'outerposition',[  -7         180        1072         909]);
    end
    saveas(gcf,figname,'fig');
%     save(fname,'VAF','Rsq'); % save for each subj
%     VAF_all(subjind,:) = nanmean(VAF); % save subj mean
%     Rsq_all(subjind,:) = nanmean(Rsq); % save subj mean
    
    close all;
end

save('force_VAF_Rsq','VAF_all','Rsq_all');
