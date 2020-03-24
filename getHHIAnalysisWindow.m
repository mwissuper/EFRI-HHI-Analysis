function [start_idx,stop_idx] = getHHIAnalysisWindow(Markers, max_distance)
%% GETHHIANALYSISWINDOW: Start and Stop indices for analyzing HHI data
%
%   [start,stop] = getHHIAnalysisWindow(Markers,max_distance) takes a
%   structure MARKERS containing VICON marker data and a maximum distance
%   MAX_DISTANCE and returns the start and stop indices for analyzing the
%   data. START and STOP are integer indices indicating the first and last
%   samples to analyze.

%   The start time is computed from the vertical displacement of the left
%   heel. Start time is the time of the first half-maximum of the heel
%   vertical displacement. Time of half-maximum is computed by finding all
%   maxima with prominence above one-fifth the range of the displacement,
%   and then finding the location of the first point greater than the first
%   half-maximum.
%   
%   Stop time is computed from the clavicle marker. The trial ends when the
%   clavicle has travelled the distance indicated in MAX_DISTANCE.

% Need to convert all to m

% Calculate the start time from the left heel's first peak vertical pos
left_heel = Markers.LHEE(:,3);
clav_Y = Markers.CLAV(:,2);
peak = findpeaks(left_heel, 'MinPeakProminence',range(left_heel)/5);
start_idx = find(left_heel >= peak(1)/2, 1, 'first');

% To find the END time, check the clavicle marker for
% displacement in accordance with MAX_DISTANCE
clav_Y = clav_Y(start_idx:end) - clav_Y(start_idx);
stop_idx = find(clav_Y > max_distance, 1, 'first');
if isempty(stop_idx)
    stop_idx = length(clav_Y);
end
stop_idx = stop_idx + start_idx - 1;
end