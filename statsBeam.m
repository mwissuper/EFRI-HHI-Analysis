% Run all stats tests. Use all data (early and late trials)

clear all; close all; clc;

subj_array = 3:14;
subj_array_force = [3:5 8:13]; % HHI_01 and HHI_02 are pilots. HHI_06, _07, _14 are all missing Fx in the force data (force was unplugged)

%% Performance metrics

kinem = load('HHI2017_Stats_MW'); % More participants than force data
conds = {'Solo Beam','Assist Beam','Solo Ground','Assist Ground'};
n = 0;
for subj = subj_array
    n = n + 1;
    for i = 1:length(conds)
        rows = kinem.GroupStats.Subject==subj & strcmp(kinem.GroupStats.Type,conds{i});
        StdSway(n,i) = nanmean(kinem.GroupStats.StdSway(rows)); % SD Clav sway
        if i == 1 || i == 2
            Dist(n,i) = nanmean(kinem.GroupStats.Dist(rows)); % Total distance traveled on beam
            AvgSpeed(n,i) = nanmean(kinem.GroupStats.AvgSpeed(rows)); % Average speed of walking on beam
        end
    end
end

%% Stats tests to compare Solo Beam vs. Assist Beam

% Compare value vs. mean (full beam length for all assisted trials)
[Hl,P1] = lillietest(Dist(:,1)); % First test if distrib is normal so know which test to use next
if P1 > 0
    [H,p1,CI,stats] = ttest(Dist(:,1),3.65);
    test1 = 't-test';
end

[p(1),stat(1),test{1}] = comp2groups(AvgSpeed(:,1),AvgSpeed(:,2));
[p(2),stat(2),test{2}] = comp2groups(StdSway(:,1),StdSway(:,2));

[stats.tstat p1]
[stat' p']
test

%% Stats tests to compare Solo Ground vs. Assist Ground
clc
% Compare value vs. mean (full beam length for all assisted trials)
[p,stat,test] = comp2groups(StdSway(:,3),StdSway(:,4));

%% Plot correlation solo balance ability (as measured by solo beam
% distance avg all trials) vs. improvement in performance (sway, beam distance, and avg
% speed) with partner. Test correlation.

kinem = load('HHI2017_Stats_MW'); % More participants than force data
conds = {'Solo Beam'};
n = 0;
for subj = subj_array
    n = n + 1;
    for i = 1:length(conds)
        rows = kinem.GroupStats.Subject==subj & strcmp(kinem.GroupStats.Type,conds{i});
        BaseDist(n,i) = nanmean(kinem.GroupStats.Dist(rows)); % Total distance traveled on beam
        a = kinem.GroupStats.Dist(rows);
        clear a
    end
end

% Calculate correlations

dSway = diff(StdSway')'; dSpeed = diff(AvgSpeed')';

% Characterizing baseline
% balance ability by mean of distance completed on beam across all late
% trials

[rho(1),p(1)] = corr(Dist(:,1),dSpeed);
[rho(2),p(2)] = corr(Dist(:,1),dSway);

%% Stats tests to compare solo ground vs. assist ground

[p,stat,test] = comp2groups(StdSway(:,3),StdSway(:,4))
[p,stat,test] = comp2groups(AvgSpeed(:,3),AvgSpeed(:,4))

%% Calc means force, power, and arm length metrics compare beam vs. ground

forces = load('HHI2017_Stats_force_MW');
conds = {'Assist Ground','Assist Beam'};
kx = []; bx = []; PposXarray = []; PnegXarray = [];
% Concatenate data
n = 0;
for subj = subj_array_force
    n = n + 1;
    for i = 1:length(conds)
        rows = forces.GroupStats.Subject==subj & strcmp(forces.GroupStats.Type,conds{i});
        
        % Force (signed)
        Fx(n,i) = nanmean(forces.GroupStats.meanFx(rows)); 
        Fy(n,i) = nanmean(forces.GroupStats.meanFy(rows));
        Fz(n,i) = nanmean(forces.GroupStats.meanFz(rows));
        FxSD(n,i) = nanmean(forces.GroupStats.SDFx(rows)); 
        FySD(n,i) = nanmean(forces.GroupStats.SDFy(rows));
        FzSD(n,i) = nanmean(forces.GroupStats.SDFz(rows));
        
        %% IP
        % Vel mag (unsigned)
        Vx(n,i) = nanmean(forces.GroupStats.meanVx(rows)); 
        SDVx(n,i) = nanmean(forces.GroupStats.SDVx(rows)); 
        Vy(n,i) = nanmean(forces.GroupStats.meanVy(rows));
        Vz(n,i) = nanmean(forces.GroupStats.meanVz(rows));
        
        % Power (signed)
        Px(n,i) = nanmean(forces.GroupStats.meanIPpowerX(rows)); 
        Py(n,i) = nanmean(forces.GroupStats.meanIPpowerY(rows));
        Pz(n,i) = nanmean(forces.GroupStats.meanIPpowerZ(rows));
        PxSD(n,i) = nanmean(forces.GroupStats.SDIPpowerX(rows)); 
        PySD(n,i) = nanmean(forces.GroupStats.SDIPpowerY(rows));
        PzSD(n,i) = nanmean(forces.GroupStats.SDIPpowerZ(rows));
        
        % effective arm length
        armPOBx(n,i) = nanmean(forces.GroupStats.meanArmPOBX(rows)); 
        SDarmPOBx(n,i) = nanmean(forces.GroupStats.SDarmPOBX(rows)); 
        
        if i == 2
            % Do only for assist beam cond. Want to compare magnitude, so
            % take abs of neg force/power
            FposX(n) = nanmean(forces.GroupStats.meanPosFx(rows)); 
            FnegX(n) = nanmean(abs(forces.GroupStats.meanNegFx(rows)));
            FposY(n) = nanmean(forces.GroupStats.meanPosFy(rows)); 
            FnegY(n) = nanmean(abs(forces.GroupStats.meanNegFy(rows)));
            FposZ(n) = nanmean(forces.GroupStats.meanPosFz(rows)); 
            FnegZ(n) = nanmean(abs(forces.GroupStats.meanNegFz(rows)));
            
            PposX(n) = nanmean(forces.GroupStats.meanPosPowerIntPtX(rows)); 
            PnegX(n) = nanmean(abs(forces.GroupStats.meanNegPowerIntPtX(rows))); 
            PposY(n) = nanmean(forces.GroupStats.meanPosPowerIntPtY(rows)); 
            PnegY(n) = nanmean(abs(forces.GroupStats.meanNegPowerIntPtY(rows)));
            PposZ(n) = nanmean(forces.GroupStats.meanPosPowerIntPtZ(rows)); 
            PnegZ(n) = nanmean(abs(forces.GroupStats.meanNegPowerIntPtZ(rows)));
            
            % Check effect of F drift on signed power
            meanPposFlo(n) = nanmean(forces.GroupStats.meanIPpowerPosFlo(rows));
            meanPnegFlo(n) = nanmean(abs(forces.GroupStats.meanIPpowerNegFlo(rows)));
            meanPposFhi(n) = nanmean(forces.GroupStats.meanIPpowerPosFhi(rows));
            meanPnegFhi(n) = nanmean(abs(forces.GroupStats.meanIPpowerNegFhi(rows)));
            
            % Percent of time doing strategy
            perPpos(n) = nanmean(forces.GroupStats.perPposFposX(rows) + forces.GroupStats.perPposFnegX(rows));
            perFpos(n) = nanmean(forces.GroupStats.perPposFposX(rows) + forces.GroupStats.perPnegFposX(rows));
            perPneg(n) = nanmean(forces.GroupStats.perPnegFposX(rows) + forces.GroupStats.perPnegFnegX(rows));
            perFneg(n) = nanmean(forces.GroupStats.perPposFnegX(rows) + forces.GroupStats.perPnegFnegX(rows)) ;
            
            % Concat stiffness and damping values, don't take mean across
            % trials! Do same for pos and neg P
            kx = [kx; forces.GroupStats.kx_torso(rows)];
            bx = [bx; forces.GroupStats.bx_torso(rows)];
            PposXarray = [PposXarray; forces.GroupStats.meanPosPowerIntPtX(rows)];
            PnegXarray = [PnegXarray; forces.GroupStats.meanNegPowerIntPtX(rows)];
            
            % Concat stiffness and damping for k-means clustering. Do mass
            % too just in case
            coeffs_array(n,:) = [nanmean(forces.GroupStats.mx_torso(rows)) nanmean(forces.GroupStats.bx_torso(rows)) nanmean(forces.GroupStats.kx_torso(rows))];
            
            % Calculate forces
            % Concat lag F to torso state
            lagFIPTorsoX(n) = nanmean(forces.GroupStats.lagFIPTorsoX(rows));
            lagFIPvTorsoX(n) = nanmean(forces.GroupStats.lagFIPvTorsoX(rows));
        end
        
        %% POB torso
        % Power mag (unsigned)
        Px_POB(n,i) = nanmean(forces.GroupStats.meanAbsPowerPOBX(rows)); 
        SDPx_POB(n,i) = nanmean(forces.GroupStats.SDAbsPowerPOBX(rows)); 
        if i == 2
            PposX_POB(n) = nanmean(forces.GroupStats.meanPosPowerPOBX(rows)); 
            PnegX_POB(n) = nanmean(abs(forces.GroupStats.meanNegPowerPOBX(rows))); 
             % Percent of time doing strategy
%             perPpos_POB(n) = nanmean(forces.GroupStats.POBperPposFposX(rows) + forces.GroupStats.POBperPposFnegX(rows));
%             perFpos_POB(n) = nanmean(forces.GroupStats.POBperPposFposX(rows) + forces.GroupStats.POBperPnegFposX(rows));
%             perPneg_POB(n) = nanmean(forces.GroupStats.POBperPnegFposX(rows) + forces.GroupStats.POBperPnegFnegX(rows));
%             perFneg_POB(n) = nanmean(forces.GroupStats.POBperPposFnegX(rows) + forces.GroupStats.POBperPnegFnegX(rows)) ;   
        end
    end
end

%% Force (mean) tests
% Compare Assist Beam vs. Light Touch
% [px1,stat,test] = compMean(Fx(:,2),1);
% [py1,stat,test] = compMean(Fy(:,2),1);
% [pz1,stat,test] = compMean(Fz(:,2),1);

% % Compare Assist Ground vs. Sylos-Labini handholding walking (they reported
% % 2-3 N, so just compare vs. 3N)
% [px,stat,test] = compMean(Fx(:,1),3)
% [py,stat,test] = compMean(Fy(:,1),3)
% [pz,stat,test] = compMean(Fz(:,1),3)
clc
clear stat 
test = {};
% Compare Assist Ground vs. Assist Beam per direction
% Stats tests compare 2 paired samples
[px,stat(1),test{1}] = comp2groups(Fx(:,1),Fx(:,2));
% [py,stat(2),test{2}] = comp2groups(Fy(:,1),Fy(:,2));
% [pz,stat(3),test{3}] = comp2groups(Fz(:,1),Fz(:,2));

[stat' [px;py;pz]]
test

clear stat 
test = {};
% For Assist Beam, compare pos vs. neg force
[px,stat(1),test{1}] = comp2groups(FposX,FnegX);
% [py,stat(2),test{2}] = comp2groups(FposY,FnegY);
% [pz,stat(3),test{3}] = comp2groups(FposZ,FnegZ);
[stat' [px;py;pz]]
test

% Worst case scenario of drift 1.5N to reduce mean differences
[px,stat,test] = comp2groups(FposX-1.5,FnegX+1.5)

%% Compare force means for better vs. worse halves (based on solo dist
% metric). 1) divide the two groups of better vs. worse based on all 12
% subj's in kinem data. 2) divide two groups based on 9 subj's with F data
% and disregard subj that was in the middle (HHI 10)
clc
% W = [9 5 3 11 10 12]; B = [13 8 4]; % Not enough elements in B to run lillietest
% for i = 1:length(W)
%     indW(i) = find(subj_array_force==W(i),1,'first');
% end
% for i = 1:length(B)
%     indB(i) = find(subj_array_force==B(i),1,'first');
% end
% [p1,stat1,test1] = comp2groups(Fx(indW,2),Fx(indB,2))
%
clear W B indW indB p stat test
[c,ia,ib] = intersect(subj_array,subj_array_force);
soloDistF = Dist(ia,1);
[soloDistFSort,iSort] = sort(soloDistF);
W = subj_array_force(iSort(1:4)); B = subj_array_force(iSort(6:end)); 

[p,stat,test] = comp2groups(Fx(iSort(1:4),2),Fx(iSort(6:end),2))

% Plot for above test
plot(1:2,[Fx(iSort(1:4),2) Fx(iSort(6:end),2)],'Fx')
xlim([0.5 2.5]),set(gca,'xtick',1:2); 
xlab{1} = 'Worse'; xlab{2} = 'Better'; set(gca,'xticklabel',xlab)
xlabel('Solo balance (beam distance)'); ylabel('Mean F mag (N)');
box off,set(gca,'tickdir','out')

%% Test correlation better/worse solo balance vs mean 
[c,ia,ib] = intersect(subj_array,subj_array_force);
[r, p] = corr(Dist(ia),Fx(:,2));

%% Force SD tests compare Assist Ground vs. Assist Beam
clear stat; test = {}; clc;
[px,stat(1),test{1}] = comp2groups(FxSD(:,1),FxSD(:,2));
[py,stat(2),test{2}] = comp2groups(FySD(:,1),FySD(:,2));
[pz,stat(3),test{3}] = comp2groups(FzSD(:,1),FzSD(:,2));

[stat' [px;py;pz]]
test

%% Calculate abs and percent change assist ground vs assist beam

deltaFx = mean(Fx(:,2))-mean(Fx(:,1))
deltaSDFx = mean(FxSD(:,2))-mean(FxSD(:,1))
% deltaFz = mean(Fz(:,2))-mean(Fz(:,1))

% deltaFxPer = (mean(Fx(:,2))-mean(Fx(:,1)))/mean(Fx(:,1))
% deltaFzPer = (mean(Fz(:,2))-mean(Fz(:,1)))/mean(Fz(:,1))

%% Velocity magnitude tests
clc
clear stat 
test = {};

% Compare Assist Ground vs. Assist Beam per direction
% Stats tests compare 2 paired samples
[px,stat(1),test{1}] = comp2groups(Vx(:,1),Vx(:,2));
[p2,stat(2),test{2}] = comp2groups(SDVx(:,1),SDVx(:,2));
% [py,stat(2),test{2}] = comp2groups(Vy(:,1),Vy(:,2));
% [pz,stat(3),test{3}] = comp2groups(Vz(:,1),Vz(:,2));

[stat' [px;p2]]
test

%% Mean power tests - IP
clc;
% Compare Assist Beam vs. normal walking (or other small-power) thresh

% clear stat 
% test = {};
% % Compare Assist Beam per direction vs. means for one person holding sensor
% % (HHI08 Assist Solo trials)
% [px,stat(1),test{1}] = compMean(Px(:,2),0.0256);
% [py,stat(2),test{2}] = compMean(Py(:,2),2.621);
% [pz,stat(3),test{3}] = compMean(Pz(:,2),0.2597);
% [stat' [px;py;pz]]
% test

clear stat 
test = {};
% Compare Assist Ground vs. Assist Beam per direction
% Stats tests compare 2 paired samples
[px,stat(1),test{1}] = comp2groups(Px(:,1),Px(:,2));
[py,stat(2),test{2}] = comp2groups(Py(:,1),Py(:,2));
[pz,stat(3),test{3}] = comp2groups(Pz(:,1),Pz(:,2));
[stat' [px;py;pz]]
test

clear stat 
test = {};

%% For Assist Beam, compare pos vs. neg power
clc; clear stat; test = {};

[px,stat(1),test{1}] = comp2groups(PposX,PnegX);
[py,stat(2),test{2}] = comp2groups(PposY,PnegY);
[pz,stat(3),test{3}] = comp2groups(PposZ,PnegZ);
[stat' [px;py;pz]]
test

% Check effect of F drift

[plo,statlo,test] = comp2groups(meanPposFlo,meanPnegFlo)

[phi,stathi,test] = comp2groups(meanPposFhi,meanPnegFhi)

%% Power SD tests compare Assist Ground vs. Assist Beam
clear stat; test = {}; clc;
[px,stat(1),test{1}] = comp2groups(PxSD(:,1),PxSD(:,2));
[py,stat(2),test{2}] = comp2groups(PySD(:,1),PySD(:,2));
[pz,stat(3),test{3}] = comp2groups(PzSD(:,1),PzSD(:,2));

[stat' [px;py;pz]]
test

%% Calculate abs and percent change assist ground vs assist beam
% deltaPx = mean(Px(:,2))-mean(Px(:,1))
deltaPy = mean(Py(:,2))-mean(Py(:,1))

% deltaPxPer = (mean(Px(:,2))-mean(Px(:,1)))/mean(Px(:,1))
% deltaPyPer = (mean(Pz(:,2))-mean(Py(:,1)))/mean(Py(:,1))

%% Calculate abs and percent difference in pos vs. negative power for assist beam
dPxBeam = mean(PnegX)-mean(PposX)

dPyBeam = mean(PnegY)-mean(PposY)

dPzBeam = mean(PnegZ)-mean(PposZ)

%% Test percent time of trials spent in each strategy

clear stat p; test = {};
[p(1),stat(1),test{1}] = comp2groups(perFpos,perFneg);
[p(2),stat(2),test{2}] = comp2groups(perPpos,perPneg);

[stat' p']
test

%% Power tests - POB torso
clc;
% Compare Assist Beam vs. normal walking 
clear stat text p
test = {};

% Compare Assist Ground vs. Assist Beam per direction
% Stats tests compare 2 paired samples
text(1,:) = 'Mean POB Px       ';
text(2,:) = 'SD POB Px         ';
text(3,:) = 'Pos Neg POB Fx    ';
text(4,:) = 'Per Pos Neg POB Px';
[p(1),stat(1),test{1}] = comp2groups(Px_POB(:,1),Px_POB(:,2));
[p(2),stat(2),test{2}] = comp2groups(SDPx_POB(:,1),SDPx_POB(:,2));
% For Assist Beam, compare pos vs. neg power
[p(3),stat(3),test{3}] = comp2groups(PposX_POB,PnegX_POB);
[p(4),stat(4),test{4}] = comp2groups(perPpos_POB,perPneg_POB);
t = table(text,test',stat',p');

%% Correlation between metrics of balance assistance and sway variability improvement - what is helping balance?
% see tests in plotGroupData with plots
%% Test effective POB arm length x dir

clear stat p; test = {};
[p(1),stat(1),test{1}] = comp2groups(armPOBx(:,1),armPOBx(:,2));
[p(2),stat(2),test{2}] = comp2groups(SDarmPOBx(:,1),SDarmPOBx(:,2));

[stat' p']
test

%% Correlation kx, bx and PosX and PnegPos using all trials b/c fits varied so much trial to trial
% remove all nan values before do corr
clear r p
indk = find(isnan(kx)==0);
indb = find(isnan(bx)==0);
[r(1),p(1)] = corr(kx(indk),PposXarray(indk));
[r(2),p(2)] = corr(kx(indk),-PnegXarray(indk));
[r(3),p(3)] = corr(bx(indb),PposXarray(indb)); % damper can only contribute to energy dissipation

[p' r']
% Plot for above corr
subplot(1,3,1),plot(kx(indk),PposXarray(indk),'kx')
xlabel('Stiffness (N/m)'),ylabel('Braking power (W)');
box off,set(gca,'tickdir','out');

subplot(1,3,2),plot(kx(indk),-PnegXarray(indk),'kx'),xlabel('Stiffness (N/m)'), hold on
c = polyfit(kx(indk),-PnegXarray(indk),1); plot(kx(indk),polyval(c,kx(indk)),'k');
ylabel('Motor power (W)');
box off,set(gca,'tickdir','out'),titlename = sprintf('rho = %.2f',r(2));
title(titlename);

subplot(1,3,3),plot(bx(indb),PposXarray(indb),'kx'),xlabel('Damping (N/(m/s))')
ylabel('Braking power (W)');
box off,set(gca,'tickdir','out')

%% Test if lags between F and torso state signif

clear stat p; test = {};
% Compare value vs. zero
[Hl,P1] = lillietest(lagFIPTorsoX); % First test if distrib is normal so know which test to use next
if P1 > 0
    [H,p,CI,stats] = ttest(lagFIPTorsoX,0);
    test = 't-test';
end

% Compare value vs. zero
[Hl,P1] = lillietest(lagFIPvTorsoX); % First test if distrib is normal so know which test to use next
if P1 > 0
    [H,p,CI,stats] = ttest(lagFIPvTorsoX,0);
    test = 't-test';
end

%% Split up data into groups based on regression coefficients, test if 
% performance improvement is sig diff between groups - one group only has 2
% partnerships in it

% Just k and b
X = coeffs_array(:,2:3);
rng(1); % For reproducibility
[idx,C] = kmeans(X,2);

% Use kmeans to compute the distance from each centroid to points on a grid. To do this, pass the centroids (C) and points on a grid to kmeans, and implement one iteration of the algorithm.
x1 = min(X(:,1)):0.01:max(X(:,1));
x2 = min(X(:,2)):0.01:max(X(:,2));
[x1G,x2G] = meshgrid(x1,x2);
XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot

idx2Region = kmeans(XGrid,2,'MaxIter',1,'Start',C);
% kmeans displays a warning stating that the algorithm did not converge, which you should expect since the software only implemented one iteration.

% Plot the cluster regions.

figure;
gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
    [0,0.75,0.75;0.75,0,0.75;0.75,0.75,0],'..');
hold on;
plot(X(:,1),X(:,2),'k*','MarkerSize',5);
xlabel('b (N/(m/s))');ylabel('k (N/m)');
hold off;

%% Split up data into groups based on R^2 each component
load('force_VAF_Rsq.mat')

% Just k and b
X = Rsq_all(:,2:3);
rng(1); % For reproducibility
[idx,C] = kmeans(X,2);

% Use kmeans to compute the distance from each centroid to points on a grid. To do this, pass the centroids (C) and points on a grid to kmeans, and implement one iteration of the algorithm.
x1 = min(X(:,1)):0.01:max(X(:,1));
x2 = min(X(:,2)):0.01:max(X(:,2));
[x1G,x2G] = meshgrid(x1,x2);
XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot

idx2Region = kmeans(XGrid,2,'MaxIter',1,'Start',C);
% kmeans displays a warning stating that the algorithm did not converge, which you should expect since the software only implemented one iteration.

% Plot the cluster regions.

figure;
gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
    [0,0.75,0.75;0.75,0,0.75;0.75,0.75,0],'..');
hold on;
plot(X(:,1),X(:,2),'k*','MarkerSize',5);
xlabel('R^2 Fb');ylabel('R^2 Fk');
hold off;

%% Test linear correlations
conds = {'Solo Beam','Assist Beam','Solo Ground','Assist Ground'};

%% Sway reduction solo to partner beam-walking
[c,ia,ib] = intersect(subj_array,subj_array_force);
SwayVRedF = StdSway(ia,1) - StdSway(ia,2);
plotind = 0;
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(Fx(:,2),SwayVRedF);
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(FxSD(:,2),SwayVRedF);
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(Vx(:,2),SwayVRedF);
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(SDVx(:,2),SwayVRedF);
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(Px(:,2),SwayVRedF);
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(PxSD(:,2),SwayVRedF);

%% More correlations 
[c,ia,ib] = intersect(subj_array,subj_array_force);
plotind = 0;
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(Fx(:,2),StdSway(ia,2));
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(Fx(:,2),coeffs_array(:,2)); % b
plotind = plotind+1;
[r(plotind),p(plotind)] = corr(Fx(:,2),coeffs_array(:,3)); % k
% plotind = plotind+1;
% [r(plotind),p(plotind)] = corr(FxSD(:,2),StdSway(ia,2));

