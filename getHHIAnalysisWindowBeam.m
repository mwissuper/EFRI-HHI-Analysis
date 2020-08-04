function [start_idx,stop_idx,beamHt] = getHHIAnalysisWindowBeam(Markers,vLHEEZfilt)
%% GETHHIANALYSISWINDOW: Start and Stop indices for analyzing HHI data
%
%   [start,stop] = getHHIAnalysisWindow(Markers,max_distance) takes a
%   structure MARKERS containing VICON marker data and returns the start and stop indices for analyzing the
%   data. START and STOP are integer indices indicating the first and last
%   samples to analyze.

%   The start time is computed from the vertical displacement of the left
%   heel. Start time is the time of the first maximum of the heel
%   vertical displacement. Maxima are computed by finding all
%   maxima with prominence above one-fifth the range of the displacement.
%   
%   Stop time is computed from vertical displacement of the left heel, when
%   both heels have dropped below beam height and stay below beam height
%   for remainder of trial

% Need to convert all to m

% Calculate the start time from the left heel's first peak vertical pos
left_heel = Markers.LHEE(:,3);
[peaks,ipk] = findpeaks(left_heel, 'MinPeakProminence',range(left_heel)/5);
start_idx = ipk(1); % find(left_heel >= peaks(1)/2, 1, 'first'); % If poor gapfill beg of data then affects finding start index. Peak should be safe though.

right_heel = Markers.RHEE(:,3);

% LHEE velocity gets close to zero between first two LHEE peaks
% vthresh = 0.01; % (mm/s)
if length(ipk) > 1
    [troughs,iTroughs] = findpeaks(-left_heel(ipk(1):ipk(2)),'MinPeakProminence',range(left_heel)/5);
else
    [troughs,iTroughs] = findpeaks(-left_heel(ipk(1):end),'MinPeakProminence',range(left_heel)/5);
end
if isempty(troughs) % sometimes the trough is too flat to be picked out by findpeaks, so look at abs vel LHEE vertical close to zero for a while
    if length(ipk) > 1
        iTroughs = find(abs(vLHEEZfilt(ipk(1):ipk(2))) < 0.05,10,'first');
    else
        iTroughs = find(abs(vLHEEZfilt(ipk(1):end)) < 0.05,10,'first');
    end
end

% To find the END time, find the first time the LHEE or RHEE hits
% ground (local min) after the first time it drops below beam height
% after the start time and stays below the beam height from that point
% onwards
indHeight = iTroughs(1)+ipk(1)-1;
beamHt = left_heel(indHeight); % for plotting ref, in mm. % Estimate beam height from height of LHEE marker at first time vertical
stop1 = find(left_heel(start_idx:end) >= beamHt, 1, 'last') + start_idx - 1;
stop2 = find(right_heel(start_idx:end) >= beamHt, 1, 'last') + start_idx - 1;

% stop as first time heel drops belows beam ht and stays there for rest
% of trial. issues with additional steps after trial ended and heel
% went up so stop is too late. Just fix manually.
stop_idx = min([stop1 stop2]);
% Below is attempt at algo to find when heel hit floor after off beam
% instead of when heel first dropped below beam ht. However, it doesn't
% work very consistently. Just go through and pick by eye. if
% min([stop1 stop2]) >= length(left_heel)-2 % can't find a peak if not
% >= 3 samples left, happens when not find beamHt correctly (need to
% fix manually in mainWorkPowerAnalysisMW code)
%     stop_idx = min([stop1 stop2]);
% else
%     if stop1 < stop2 % look for min on left side
%         [troughs,iTroughs] =
%         findpeaks(-left_heel(stop1:end),'MinPeakProminence',range(left_heel)/50);
%     else
%         [troughs,iTroughs] =
%         findpeaks(-right_heel(stop2:end),'MinPeakProminence',range(left_heel)/50);
%     end if isempty(iTroughs) % too flat, no obvious troughs for HS
%     event
%         stop_idx = min([stop1 stop2]);
%     else
%         stop_idx = iTroughs(1) + min([stop1 stop2]) - 1;
%     end
% end

% % Plots to check plot(left_heel) hold on; plot(right_heel)
% vline([start_idx stop_idx],'k-') vline(indHeight,'g--','on beam');
% hline(beamHt,'k-'); plot(ipk,peaks,'x') plot(iTroughs(1)+min([stop1
% stop2])-1,-troughs,'go') % plot when foot stepped onto floor
% ylabel('Vertical pos (mm)');






