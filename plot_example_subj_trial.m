% Find example subject using LateStats mat file
% First run plotGroupData.m to after line 120

i = 2;

rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds{i});

%% Sway - find subj with Solo and Partner beam trials closest to group mean
conds = {'Solo Beam','Assist Beam'};
[m,imin] = min(abs(StdSway-nanmean(StdSway)));

subj_array_force(imin);

% After choose which subj to use
subj = 8; % Subj 8 is close to mean for ML and Vert AG and AB cond's

for i = 3%:4 % Solo beam and partner beam cond's
    clear imin row trial
    rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds{i-2});
    SwaySubj = LateGroup.StdSway(rows);
    [m,imin] = min(abs(SwaySubj-nanmean(StdSway(:,i)))); % Find indiv trial with value closes to group mean
    ind = find(rows==1,imin,'first');
    row{i-2} = ind(imin);
    trial(i-2) = LateGroup.TrialNumber(row{i-2});
end

% HHI08 T11 mean sway close to group mean sway for partner beam, T42 close
% to group mean for solo beam

%% Fz
conds = {'Assist Ground','Assist Beam'};
[m,imin] = min(abs(meanFz-nanmean(meanFz)));

subj_array_force(imin)

% After choose which subj to use
subj = 8; % Subj 8 is close to mean for ML and Vert AG and AB cond's
rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds{i});
FzSubj = LateGroup.meanFz(rows);
[m,imin] = min(abs(FzSubj-nanmean(meanFz(:,i))));

ind = find(rows==1,imin,'first');

row = ind(imin);
trial = LateGroup.TrialNumber(row);

% Fz examples: HHI8, trial 32 (AG) and trial 47 (AB)