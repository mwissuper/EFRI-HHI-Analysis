function [start_idx,stop_idx,ht,midline] = getHHIAnalysisWindow_MW(Markers,vLHEEZfilt,beam,torso,subj,trial)
%% GETHHIANALYSISWINDOW: Start and Stop indices for analyzing HHI data
%
%   [start,stop] = getHHIAnalysisWindow(Markers,max_distance) takes a
%   structure MARKERS containing VICON marker data and returns the start and stop indices for analyzing the
%   data. START and STOP are integer indices indicating the first and last
%   samples to analyze.

%   The start time is computed from the vertical displacement and vel of the left
%   heel. Start time is the time of the first maximum of the heel
%   vertical displacement. Maxima are computed by finding all
%   maxima with prominence above one-fifth the range of the displacement.
%   
%   Stop time is computed from vertical displacement of the left heel, when
%   both heels have dropped below beam height and stay below beam height
%   for remainder of trial. When no beam, look at torso forward vel close
%   to zero.

beamLen = 3.4; % (m), subtract off some distance from actual length of 3.65m so not capture lots of standing time at end of trial

% Calculate the start time from the left heel's first peak vertical pos
left_heel = Markers.LHEE/1000;
right_heel = Markers.RHEE/1000;
right_toe = Markers.RTOE/1000;

if beam == 1 % Sometimes small peaks so use different thresh than no beam
    [peaks,ipk] = findpeaks(left_heel(:,3), 'MinPeakProminence',range(left_heel(:,3))/5);
    start_idx = ipk(1); 

    % LHEE velocity gets close to zero between first two LHEE peaks
    % vthresh = 0.01; % (mm/s)
    if length(ipk) > 1
        [troughs,iTroughs] = findpeaks(-left_heel(ipk(1):ipk(2),3),'MinPeakProminence',range(left_heel(:,3))/5);
    else
        [troughs,iTroughs] = findpeaks(-left_heel(ipk(1):end,3),'MinPeakProminence',range(left_heel(:,3))/5);
    end
    if isempty(troughs) % sometimes the trough is too flat to be picked out by findpeaks, so look at abs vel LHEE vertical close to zero for a while
        if length(ipk) > 1
            iTroughs = find(abs(vLHEEZfilt(ipk(1):ipk(2))) < 0.05,10,'first');
        else
            iTroughs = find(abs(vLHEEZfilt(ipk(1):end)) < 0.05,10,'first');
        end
    end
    % To find the END time, find the first time the LHEE or RHEE hits ground (local min) after the first time it drops below beam
    % height after the start time and stays below the beam height from that
    % point onwards
    
    indHeight = iTroughs(1)+ipk(1)-1; % Index at which to get beam height
    
    % special exceptions to beamHt due to errors in gapfill preprocessing
    if subj == 3 
        if trial == 2 || trial ==4 
            ht = left_heel(464,3);
        elseif trial == 27 || trial == 31
            ht = left_heel(600,3);
        elseif trial == 7 || trial == 37 || trial == 44 % never got left foot on beam
            ht = right_heel(indHeight,3);
        else
            ht = left_heel(indHeight,3); % for plotting ref, in mm. % Estimate beam height from height of LHEE marker at first time vertical
        end
    elseif subj == 4 && trial == 1
        ht = left_heel(635,3);
    elseif subj == 5
        if trial == 3 || trial == 7
            ht = left_heel(600,3);
        else
            ht = left_heel(indHeight,3);
        end
     elseif subj == 7
        if trial == 7
            ht = left_heel(400,3);
        elseif trial == 20
            ht = right_heel(250,3);
        elseif trial == 21
            ht = left_heel(330,3);
        else
            ht = left_heel(indHeight,3);
        end
    elseif subj == 9
        if trial == 1
            ht = left_heel(400,3);
        elseif trial == 4 
            ht = right_heel(200,3);
        elseif trial == 5
            ht = right_heel(300,3);
        else
            ht = left_heel(indHeight,3);
        end
    elseif subj == 12
        if trial == 3 || trial == 26
            ht = right_heel(200,3);
        elseif trial == 43
            ht = right_toe(250,3);
        else
            ht = left_heel(indHeight,3);
        end
    elseif subj == 13
        if trial == 9 || trial == 34
            ht = right_heel(200,3);
        elseif trial == 17
            ht = right_toe(250,3);
        else
            ht = left_heel(indHeight,3);
        end
    else
        ht = left_heel(indHeight,3); % Estimate beam height from height of LHEE marker at first trough. Can't use RHEE due to gapfill errors
    end
    % stop as first time heel drops belows beam ht and stays there for rest of
    % trial. 
    stop1 = find(left_heel(start_idx:end,3) >= ht, 1, 'last') + start_idx - 1;
    stop2 = find(right_heel(start_idx:end,3) >= ht, 1, 'last') + start_idx - 1;
    stop_idx = min([stop1 stop2]);
    
    % special exceptions to midline due to marker poor fill
    if subj == 5
        if trial == 18 || trial == 49 % don't keep the value for midline
            midline = left_heel(indHeight,1);
        else
            midline = right_toe(start_idx,1);
        end
    else
        midline = right_toe(start_idx,1);
    end
else % no beam, find big peaks
    [peaks,ipk] = findpeaks(left_heel(:,3), 'MinPeakProminence',range(left_heel(:,3))/2);
    [Rpeaks,ipkR] = findpeaks(right_heel(:,3), 'MinPeakProminence',range(right_heel(:,3))/2);
    iSorted = sort([ipk; ipkR]); % sort all pks in chrono order
    start_idx = iSorted(1); % earliest event
    if (subj == 5 && (trial == 24 || trial == 40 || trial == 47 || trial == 50)) || (subj == 7 && (trial == 22 || trial == 48))% special exception took few steps, solo
        stop_idx = iSorted(end);
    else
        stop_idx = iSorted(6); % Guarantees we look at same num steps per participant 
    end
%     stop_idx = find((torsoY(start_idx:end) - torsoY(start_idx))> beamLen, 1, 'first') + start_idx + 1; 
%     if isempty(stop_idx) % use max forward pos of torso after start
%         [m,i] = max(torsoY(start_idx:end));
%         stop_idx = i + start_idx - 1;
%     end
%     vTorsoY_offset = [vTorsoY(2:end); nan];
%     stop_idx = find(vTorsoY_offset(start_idx:end) < 0.2 & vTorsoY(start_idx:end) > 0.2,1,'last') + start_idx - 1; % Choose arbitrary threshold for AP velocity of torso based on visual observation
%     if isempty(stop_idx)
%         stop_idx = length(vTorsoY)+1;
%     end
    
    ht = 0; % Assume height of zero for floor (depends on Vicon calibration)
    midline = nanmean(torso(start_idx:stop_idx,1));
end

%% Special exception trials for start/stop time window that can't be fixed algorithmically
% found index by eye in Nexus
if subj == 3 
    if trial == 2
        stop_idx = 1357;
    elseif trial == 3
        stop_idx = 419; % took 2 tries in this trial, just count first one
    elseif trial == 4
        stop_idx = 1447;
    elseif trial == 13 
        stop_idx = 441;
    elseif trial == 20
        stop_idx = 829;
    elseif trial == 29
        stop_idx = 1329;
    elseif trial == 36
        stop_idx = 1143;
    elseif trial == 37
        stop_idx = 442;
    elseif trial == 42
        stop_idx = 1376;
    elseif trial == 46
        stop_idx = 1204;
    end
elseif subj == 4
    if trial == 1
        stop_idx = 1350;
    elseif trial == 3
        stop_idx = 3792;
        start_idx = start_idx + 50; % transient at beginning due to incorrect gapfill
    elseif trial == 7
        stop_idx = 2972;
    elseif trial == 9
        stop_idx = 1773;
    elseif trial == 11
        stop_idx = 3908;
    elseif trial == 17 
        stop_idx = 1960;
    elseif trial == 22
        stop_idx = 2224;
    elseif trial == 24
        stop_idx = 2014;
    elseif trial == 26
        stop_idx = 1776;
    elseif trial == 35
        stop_idx = 1321;
    elseif trial == 39
        stop_idx = 1413;
    elseif trial == 40
        stop_idx = 1388;
    elseif trial == 43
        stop_idx = 1520;
    elseif trial == 45
        stop_idx = 1105;
    elseif trial == 46
        stop_idx = 1261;
    end
elseif subj == 5
    if trial == 3
        stop_idx = 779;
    elseif trial == 4
        stop_idx = 1537;
    elseif trial == 7
        stop_idx = 885;
    elseif trial == 13
        stop_idx = 685;
    elseif trial == 17
        stop_idx = 1179;
    elseif trial == 18
        stop_idx = 960;
    elseif trial == 21
        stop_idx = 467;
    elseif trial == 23
        stop_idx = 985;
    elseif trial == 29
        stop_idx = 824;
    elseif trial == 35
        stop_idx = 815;
    elseif trial == 36
        stop_idx = 796;
    elseif trial == 37
        stop_idx = 943;
    elseif trial == 42
        stop_idx = 534;
    elseif trial == 44
        stop_idx = 489;
    elseif trial == 46
        stop_idx = 842;
    elseif trial == 49
        stop_idx = 661;
    end
elseif subj == 7
    if trial == 3
        stop_idx = 1436;
    elseif trial == 7
        stop_idx = 1176;
    elseif trial == 20
        stop_idx = 1139;
    elseif trial == 21
        stop_idx = 997;
    elseif trial == 27
        stop_idx = 1146;
    elseif trial == 37
        stop_idx = 963;
    end
elseif subj == 8
    if trial == 5
        stop_idx = 1042;
    elseif trial == 8
        stop_idx = 987;
    elseif trial == 9
        stop_idx = 807;
    elseif trial == 10
        stop_idx = 712;
    elseif trial == 11
        stop_idx = 953;
    elseif trial == 13
        stop_idx = 807;
    elseif trial == 15
        stop_idx = 810;
    elseif trial == 21
        stop_idx = 774;
    elseif trial == 23
        stop_idx = 838;
    elseif trial == 24
        stop_idx = 764;
    elseif trial == 28
        stop_idx = 772;
    elseif trial == 36
        stop_idx = 822;
    elseif trial == 39
        stop_idx = 592;
    elseif trial == 42
        stop_idx = 634;
    elseif trial == 46
        stop_idx = 757;
    elseif trial == 47
        stop_idx = 814;
    elseif trial == 48
        stop_idx = 693;
    end
elseif subj == 9
    if trial == 1
        stop_idx = 941;
    elseif trial == 5
        stop_idx = 1659;
    elseif trial == 6
        stop_idx = 1564;
    elseif trial == 17
        stop_idx = 1265;
    elseif trial == 20
        stop_idx = 1313;
    elseif trial == 22
        stop_idx = 1396;
    elseif trial == 28
        stop_idx = 1258;
    elseif trial == 33
        stop_idx = 1173;
    elseif trial == 39
        stop_idx = 1202;
    elseif trial == 46
        stop_idx = 1212;
    elseif trial == 49
        stop_idx = 879;
    elseif trial == 50
        stop_idx = 1093;
    end
elseif subj == 10
    if trial == 3
        stop_idx = 1554;
    elseif trial == 11
        stop_idx = 807;
    elseif trial == 23
        stop_idx = 1116;
    elseif trial == 24
        stop_idx = 948;
    elseif trial == 45
        stop_idx = 1028;
    elseif trial == 47
        stop_idx = 885;
    end
elseif subj == 13
    if trial == 9
        stop_idx = 549;
    elseif trial == 10
        stop_idx = 1322;
    elseif trial == 34
        stop_idx = 558;
    end
elseif subj == 14
    if trial == 8
        stop_idx = 437;
    end
end

%% Plots to check
% close all;
% plot(left_heel(:,3))
% hold on;
% plot(right_heel(:,3))
% vline([start_idx stop_idx],'k-');
% % hline(beamHt,'k-');
% plot(ipk,peaks,'bx',ipkR,Rpeaks,'rx')
% %     plot(iTroughs(1)+min([stop1 stop2])-1,-troughs(1),'go') % plot when foot stepped onto floor
% ylabel('Vertical pos (m)');
