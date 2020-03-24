% Run all stats tests

clear all; close all; clc;

subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)

%% Performance metrics

kinem = load('HHI2017_LateStats_MW'); % More participants than force data
conds = {'Solo Beam','Assist Beam','Solo Ground','Assist Ground'};
n = 0;
for subj = subj_array
    n = n + 1;
    for i = 1:length(conds)
        rows = kinem.LateGroup.Subject==subj & strcmp(kinem.LateGroup.Type,conds{i});
        StdSway(n,i) = nanmean(kinem.LateGroup.StdSway(rows)); % SD Clav sway
        Dist(n,i) = nanmean(kinem.LateGroup.Dist(rows)); % Total distance traveled on beam
        AvgSpeed(n,i) = nanmean(kinem.LateGroup.AvgSpeed(rows)); % Average speed of walking on beam
    end
end

%% Stats tests to compare Solo Beam vs. Assist Beam

[p,stat,test] = comp2groups(StdSway(:,1),StdSway(:,2))
[p,stat,test] = comp2groups(AvgSpeed(:,1),AvgSpeed(:,2))

% Compare value vs. mean (full beam length for all assisted trials)
[Hl,P1] = lillietest(Dist(:,1)); % First test if distrib is normal so know which test to use next
if P1 > 0
    [H,P0,CI,STATS] = ttest(Dist(:,1),Dist(1,2))
end

%% Plot correlation baseline balance ability (as measured by solo beam
% distance) vs. improvement in performance (sway, beam distance, and avg
% speed) with partner. Test correlation.

dSway = diff(StdSway')'; dDist = diff(Dist')'; dSpeed = diff(AvgSpeed')';

subplot(1,3,1),plot(Dist(:,1),dSway,'x')
[rho,p] = corr(Dist(:,1),dSway);
if p < 0.05
    c = polyfit(Dist(:,1),dSway,1);
    hold on; plot(Dist(:,1),polyval(c,Dist(:,1)),'k');
    titlename = sprintf('rho = %.2f, p = %.2f',rho,p);
else
    titlename = sprintf('p = %.2f',p);
end
xlabel('Solo Beam distance (m)')
ylabel('Change in sway (m)');
title(titlename);

subplot(1,3,2),plot(Dist(:,1),dDist,'x')
[rho,p] = corr(Dist(:,1),dDist);
if p < 0.05
    c = polyfit(Dist(:,1),dDist,1);
    hold on; plot(Dist(:,1),polyval(c,Dist(:,1)),'k');
    titlename = sprintf('rho = %.2f, p = %.2f',rho,p);
else
    titlename = sprintf('p = %.2f',p);
end
xlabel('Solo Beam distance (m)')
ylabel('Change in beam distance (m)');
title(titlename);

subplot(1,3,3),plot(Dist(:,1),dSpeed,'x')
[rho,p] = corr(Dist(:,1),dSpeed);
if p < 0.05
    c = polyfit(Dist(:,1),dSpeed,1);
    hold on; plot(Dist(:,1),polyval(c,Dist(:,1)),'k');
    titlename = sprintf('rho = %.2f, p = %.2f',rho,p);
else
    titlename = sprintf('p = %.2f',p);
end
xlabel('Solo Beam distance (m)')
ylabel('Change in avg speed (m/s)');
title(titlename);

%% Stats tests to compare solo ground vs. assist ground

[p,stat,test] = comp2groups(StdSway(:,3),StdSway(:,4))
[p,stat,test] = comp2groups(AvgSpeed(:,3),AvgSpeed(:,4))

%% Force and power compare beam vs. ground

forces = load('HHI2017_LateStats_force_MW');
conds = {'Assist Ground','Assist Beam'};

% Concatenate data
n = 0;
for subj = subj_array_force
    n = n + 1;
    for i = 1:length(conds)
        rows = forces.LateGroup.Subject==subj & strcmp(forces.LateGroup.Type,conds{i});
        
        % Force mag (unsigned)
        Fx(n,i) = nanmean(forces.LateGroup.meanFx(rows)); 
        Fy(n,i) = nanmean(forces.LateGroup.meanFy(rows));
        Fz(n,i) = nanmean(forces.LateGroup.meanFz(rows));
        
        %% IP
        % Vel mag (unsigned)
        Vx(n,i) = nanmean(forces.LateGroup.meanVx(rows)); 
        Vy(n,i) = nanmean(forces.LateGroup.meanVy(rows));
        Vz(n,i) = nanmean(forces.LateGroup.meanVz(rows));
        
        % Power mag (unsigned)
        Px(n,i) = nanmean(forces.LateGroup.meanAbsPowerIntPtX(rows)); 
        Py(n,i) = nanmean(forces.LateGroup.meanAbsPowerIntPtY(rows));
        Pz(n,i) = nanmean(forces.LateGroup.meanAbsPowerIntPtZ(rows));
        if i == 2
            % Do only for assist beam cond. Want to compare magnitude, so
            % take abs of neg power
            PposX(n) = nanmean(forces.LateGroup.meanPosPowerIntPtX(rows)); 
            PnegX(n) = nanmean(abs(forces.LateGroup.meanNegPowerIntPtX(rows))); 
            PposY(n) = nanmean(forces.LateGroup.meanPosPowerIntPtY(rows)); 
            PnegY(n) = nanmean(abs(forces.LateGroup.meanNegPowerIntPtY(rows)));
            PposZ(n) = nanmean(forces.LateGroup.meanPosPowerIntPtZ(rows)); 
            PnegZ(n) = nanmean(abs(forces.LateGroup.meanNegPowerIntPtZ(rows)));
        end
        
        %% POB
        % Vel mag (unsigned)
        VxPOB(n,i) = nanmean(forces.LateGroup.meanVx_clav(rows)); 
        VyPOB(n,i) = nanmean(forces.LateGroup.meanVy_clav(rows));
        VzPOB(n,i) = nanmean(forces.LateGroup.meanVz_clav(rows));
        
        % Power mag (unsigned)
        PxPOB(n,i) = nanmean(forces.LateGroup.meanAbsPowerPOBX(rows)); 
        PyPOB(n,i) = nanmean(forces.LateGroup.meanAbsPowerPOBY(rows));
        PzPOB(n,i) = nanmean(forces.LateGroup.meanAbsPowerPOBZ(rows));
        if i == 2
            % Do only for assist beam cond. Want to compare magnitude, so
            % take abs of neg power
            PposXPOB(n) = nanmean(forces.LateGroup.meanPosPowerPOBX(rows)); 
            PnegXPOB(n) = nanmean(abs(forces.LateGroup.meanNegPowerPOBX(rows))); 
            PposYPOB(n) = nanmean(forces.LateGroup.meanPosPowerPOBY(rows)); 
            PnegYPOB(n) = nanmean(abs(forces.LateGroup.meanNegPowerPOBY(rows)));
            PposZPOB(n) = nanmean(forces.LateGroup.meanPosPowerPOBZ(rows)); 
            PnegZPOB(n) = nanmean(abs(forces.LateGroup.meanNegPowerPOBZ(rows)));
        end
    end
end

%% Calculate percent change assist ground vs assist beam
deltaFx = (mean(Fx(:,2))-mean(Fx(:,1)))/mean(Fx(:,1))
deltaFz = (mean(Fz(:,2))-mean(Fz(:,1)))/mean(Fz(:,1))

deltaPx = (mean(Px(:,2))-mean(Px(:,1)))/mean(Px(:,1))
deltaPz = (mean(Pz(:,2))-mean(Pz(:,1)))/mean(Pz(:,1))

%% Calculate percent difference in pos vs. negative power for assist beam
dPxBeam = (mean(PnegX)-mean(PposX))/mean(PposX)

%% Force tests
% Compare Assist Beam vs. Light Touch
[px1,stat,test] = compMean(Fx(:,2),1)
[py1,stat,test] = compMean(Fy(:,2),1)
[pz1,stat,test] = compMean(Fz(:,2),1)

% % Compare Assist Ground vs. Sylos-Labini handholding walking (they reported
% % 2-3 N, so just compare vs. 3N)
% [px,stat,test] = compMean(Fx(:,1),3)
% [py,stat,test] = compMean(Fy(:,1),3)
% [pz,stat,test] = compMean(Fz(:,1),3)

% Compare Assist Ground vs. Assist Beam per direction
% Stats tests compare 2 paired samples
[px,stat,test] = comp2groups(Fx(:,1),Fx(:,2))
[py,stat,test] = comp2groups(Fy(:,1),Fy(:,2))
[pz,stat,test] = comp2groups(Fz(:,1),Fz(:,2))

%% Velocity magnitude tests
% Compare Assist Ground vs. Assist Beam per direction
% Stats tests compare 2 paired samples
[px,stat,test] = comp2groups(Vx(:,1),Vx(:,2))
[py,stat,test] = comp2groups(Vy(:,1),Vy(:,2))
[pz,stat,test] = comp2groups(Vz(:,1),Vz(:,2))

disp('POB')
% POB
[px,stat,test] = comp2groups(VxPOB(:,1),VxPOB(:,2))
[py,stat,test] = comp2groups(VyPOB(:,1),VyPOB(:,2))
[pz,stat,test] = comp2groups(VzPOB(:,1),VzPOB(:,2))

%% Power tests

% Compare Assist Beam vs. normal walking (or other small-power) thresh

% Compare Assist Beam per direction vs. means for one person holding sensor
[px,stat,test] = compMean(Px(:,2),0.0256)
[py,stat,test] = compMean(Py(:,2),2.621)
[pz,stat,test] = compMean(Pz(:,2),0.2597)

% Compare Assist Ground vs. Assist Beam per direction
% Stats tests compare 2 paired samples
[px,stat,test] = comp2groups(Px(:,1),Px(:,2))
[py,stat,test] = comp2groups(Py(:,1),Py(:,2))
[pz,stat,test] = comp2groups(Pz(:,1),Pz(:,2))

% For Assist Beam, compare pos vs. neg power
[px,stat,test] = comp2groups(PposX,PnegX)
[py,stat,test] = comp2groups(PposY,PnegY)
[pz,stat,test] = comp2groups(PposZ,PnegZ)

disp('POB')
% Compare Assist Ground vs. Assist Beam per direction
% Stats tests compare 2 paired samples
[px,stat,test] = comp2groups(PxPOB(:,1),PxPOB(:,2))
[py,stat,test] = comp2groups(PyPOB(:,1),PyPOB(:,2))
[pz,stat,test] = comp2groups(PzPOB(:,1),PzPOB(:,2))

% For Assist Beam, compare pos vs. neg power
[px,stat,test] = comp2groups(PposXPOB,PnegXPOB)
[py,stat,test] = comp2groups(PposYPOB,PnegYPOB)
[pz,stat,test] = comp2groups(PposZPOB,PnegZPOB)

