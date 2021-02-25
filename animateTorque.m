%% Animate one trial of data to show orientation of torque on POB simultaneously with angular position

clear; clc; close all;

subj = 3; 
n = 34; % trial index

% Get mass of subj for ang momentum metric. Since compare solo to
% assist beam, only care subj's with dynamics data
if subj == 3
    mass = 49.3; % kg
    height = 1.60; % m
elseif subj == 4
    mass = 56.9; % kg
    height = 1.61; % m
elseif subj == 5
    mass = 73.8; % kg
    height = 1.64; % m
elseif subj == 6
    mass = 78.2; % kg
    height = 1.73; % m
elseif subj == 7
    mass = 76.0; % kg
    height = 1.82; % m
elseif subj == 8
    mass = 67.5; % kg
    height = 1.69; % m
elseif subj == 9
    mass = 75.2; % kg
    height = 1.68; % m
elseif subj == 10
    mass = 55.8; % kg
    height = 1.61; % m
elseif subj == 11
    mass = 62.0; % kg
    height = 1.64; % m
elseif subj == 12
    mass = 67.2; % kg
    height = 1.72; % m
elseif subj == 13
    mass = 54.8; % kg
    height = 1.61; % m
elseif subj == 14
    mass = 69.3; % kg
    height = 1.78; % m
end

filename = sprintf('HHI2017_%i.mat',subj);
load(filename);
plotind = 0;
sample_rate = TrialData(1).Markers.samplerate;

lFin = TrialData(n).Results.IP;
Force = TrialData(n).Results.Forces;
rGd = lFin - [TrialData(n).Results.midline 0 0];

%%
figure, hold on, xlim([-2 0.5]),ylim([-0.5 2]),axis square
for ind = 1:length(TrialData(n).Results.torso)
    hline(0,'k-'),vline(TrialData(n).Results.midline,'k-') % reference lines for visualization
    a = line([TrialData(n).Results.midline TrialData(n).Results.torso(ind,1)],[0 TrialData(n).Results.torso(ind,3)]); % plot torso to check angle   
    b = quiver(lFin(ind,1),lFin(ind,3),Force(ind,1),Force(ind,3),1/100); % force vector located at POB lFin    
    set(b,'color','r');
    theta
    c = quiver(TrialData(n).Results.midline,0,dx,dz); % moment arm vector
    set(c,'color','b');
    titlename = sprintf('HHI%i %s: Torque about beam = %.2f',subj,TrialData(n).Info.Trial,TrialData(n).Results.TyGd(ind)); title(titlename);
    legend('Beam to torso','Force');
    % line([TrialData(n).Results.midline rGd(ind,1)],[0 rGd(ind,3)]) % rGd from midline on beam/ground to lFin
    drawnow
    delete(a); delete(b)
end



