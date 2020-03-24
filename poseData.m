function Pose = poseData(Markers,idx)
%POSEDATA: Creates a matrix of connections between the marker points to redraw the
%stick figure
%
%   Pose = postData(Markers,idx) returns a data matrix Pose which, when
%   plotted as plot3(Pose(:,1),Pose(:,2),Pose(:,3)) draws a ball-and-stick
%   pose reconstruction of the subject whose marker data is stored in
%   Markers at sample index idx

%   Luke Drnach
%   Georgia Institue of Technology
%   February 26, 2017
    
%a row of NaNs causes the line between two points to NOT be drawn
skip = nan(1,3);

% NEED TO CHECK THAT THE MARKERS EXIST FIRST.

%Pull out the head markers
LFHD = Markers.LFHD(idx,:);
RFHD = Markers.RFHD(idx,:);
RBHD = Markers.RBHD(idx,:);
LBHD = Markers.LBHD(idx,:);
%Order the points to create the Head connections
HeadData = [LFHD;RFHD;RBHD;LFHD;LBHD;RBHD;RFHD;LBHD];
%NOTE: Matlab's Plot3 function will draw lines between successive points.
%For example, the sequence LFHD;RFHD;RBHD;LFHD will create a triangle whose
%vertices are LFHD, RFHD, and RBHD. This scheme is used throughout to
%recreate the connections between the marker points.

%Pull out the upper torso points
C7 = Markers.C7(idx,:);
CLAV = Markers.CLAV(idx,:);
LSHO = Markers.LSHO(idx,:);
RSHO = Markers.RSHO(idx,:);
RBAK = Markers.RBAK(idx,:);
%Make the torso connections
TorsoData = [LSHO; CLAV; C7; RBAK; skip; CLAV; RSHO];

%Draw the arms
if all(isfield(Markers,{'RELB','RFRM','RFIN','RUPA'}))
   RELB = Markers.RELB(idx,:);
   RFRM = Markers.RFRM(idx,:);
   RFIN = Markers.RFIN(idx,:);
   RUPA = Markers.RUPA(idx,:);
   RArmData = [RSHO; RELB; skip; RFRM; skip; RUPA; skip; RFIN];
else
    RArmData = skip;
end
if all(isfield(Markers,{'LELB','LFRM','LFIN','LUPA'}))
   LELB = Markers.LELB(idx,:);
   LFRM = Markers.LFRM(idx,:);
   LFIN = Markers.LFIN(idx,:);
   LUPA = Markers.LUPA(idx,:);
   LArmData = [LSHO; LELB; skip; LFRM; skip; LUPA; skip; LFIN];
else
    LArmData = skip;
end

%Pull out the hip markers
LASI = Markers.LASI(idx,:);
RASI = Markers.RASI(idx,:);
LPSI = Markers.LPSI(idx,:);
RPSI = Markers.RPSI(idx,:);
%Make the hip connections
HipData = [LASI;RASI;LPSI;RPSI;LASI];  

%Left leg points
LTHI = Markers.LTHI(idx,:);
LKNE = Markers.LKNE(idx,:);
LTIB = Markers.LTIB(idx,:);
LANK = Markers.LANK(idx,:);
LHEE = Markers.LHEE(idx,:);
LTOE = Markers.LTOE(idx,:);
%Left leg connections
LLegData = [LTHI;LKNE;LASI;LPSI;LKNE;LTIB;LANK;LHEE;LTOE;LANK;LKNE];

%Right leg points
RTHI = Markers.RTHI(idx,:);
RKNE = Markers.RKNE(idx,:);
RTIB = Markers.RTIB(idx,:);
RANK = Markers.RANK(idx,:);
RHEE = Markers.RHEE(idx,:);
RTOE = Markers.RTOE(idx,:);
%Right leg connections
RLegData = [RASI;RTHI;RKNE;RPSI;RASI;RKNE;RTIB;RANK;RHEE;RTOE;RANK;RKNE];

%Separate the individual body parts with skips. Note the Hip and Left leg
%do not have a skip: this is to complete the drawing of the hip using the
%left leg connections.
Pose = [HeadData;skip;TorsoData;skip;RArmData;skip;LArmData;skip;HipData;LLegData;skip;RLegData];

Pose(abs(Pose) < 1) = NaN;

% IF there is no output argument, make the plot
if nargout == 0
   figure()
   plot3(Pose(:,1),Pose(:,2),Pose(:,3),'Marker','o','LineWidth',1.5,'LineStyle','-','MarkerFaceColor','auto')
   set(get(gca,'Children'),'MarkerFaceColor',get(get(gca,'Children'),'Color'))
   
   axis equal
end

end

