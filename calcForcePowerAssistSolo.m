% Calculate mean power for Assist Solo trials where person walks while
% holding F/T sensor to get baseline power for no assist case. Also plot
% force vs. time and find mean force magnitude and mass-spring-damper fit

% Verified from videos that HHI8 is holding sensor during Assist Solo
% trials so use these. Out of other subjects, only have videos only for 
% HHI02 and HHI14, but they weren't usable for force analysis

% Also get mean raw FT measurements to get sense of bias

clear; clc; close all;
load HHI2017_8;

plot_power = 0;
plot_force = 0;
plot_force_raw = 1; % Just the voltage values to use for sensor bias
plot_force_fit = 0; % plot fit param's

trials = [2 7 12 17 22 25 33 34 41 49]; % Assist solo
% trials = [5 8 11 15 23 28 36 37 46 47]; % Assist beam
numcols = 5;

for dir = 1:3 % all components of force/power
    figure;
    plotind = 0;
    i = 0;
    for n = trials
        plotind = plotind + 1;
        % Plot pos/neg power, abs power, and mean abs power to check that the
        % mean abs power is a good metric
        subplot(2,5,plotind), hold on
        if plot_power == 1
            plot(TrialData(n).Results.time(2:end),TrialData(n).Results.IntPower(:,dir))
            plot(TrialData(n).Results.time(2:end),abs(TrialData(n).Results.IntPower(:,dir)))
            hline(mean(abs(TrialData(n).Results.IntPower(:,dir))),'k');
            if plotind == 1 || plotind == numcols + 1
                if dir == 1
                    ylabel('ML Power (W)')
                elseif dir == 2
                    ylabel('AP Power (W)')
                else
                    ylabel('Vert Power (W)')
                end
            end
        elseif plot_force == 1
            plot(TrialData(n).Results.time,TrialData(n).Results.Forces(:,dir))
%             plot(TrialData(n).Results.time,abs(TrialData(n).Results.Forces(:,dir)))
            hline(mean(TrialData(n).Results.Forces(:,dir)),'k');
            if plotind == 1 || plotind == numcols + 1
                if dir == 1
                    ylabel('ML Force (N)')
                elseif dir == 2
                    ylabel('AP Force (N)')
                else
                    ylabel('Vert Force (N)')
                end
            end
        elseif plot_force_raw == 1
            plot(TrialData(n).Results.time,TrialData(n).Results.rawF(:,dir))
            hline(mean(TrialData(n).Results.rawF(:,dir)),'k');
            if plotind == 1 || plotind == numcols + 1
                if dir == 1
                    ylabel('ML Force (N)')
                elseif dir == 2
                    ylabel('AP Force (N)')
                else
                    ylabel('Vert Force (N)')
                end
            end
        end
        
        % Labels
        xlabel('Time (s)');
        if plotind == 1
%             legend('signed','abs');
        end
        titlename = sprintf('Trial %i',n); title(titlename);

        % Calculate metrics and concatenate
        i = i + 1;
        meanPower(i,dir) = mean(abs(TrialData(n).Results.IntPower(:,dir)));
        meanForces(i,dir) = mean(abs(TrialData(n).Results.Forces(:,dir)));
        meanRawF(i,dir) = mean(TrialData(n).Results.rawF(:,dir));
    end
end