% Find example subject using LateStats mat file
% First run plotGroupData.m to after line 120

conds = {'Assist Ground','Assist Beam'};

i = 2;

rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds{i});

[m,imin] = min(abs(meanFz-nanmean(meanFz)));

subj_array_force(imin)

%% After choose which subj to use
subj = 8; % Subj 8 is close to mean for ML and Vert AG and AB cond's
rows = LateGroup.Subject==subj & strcmp(LateGroup.Type,conds{i});
FzSubj = LateGroup.meanFz(rows);
[m,imin] = min(abs(FzSubj-nanmean(meanFz(:,i))));

ind = find(rows==1,imin,'first');

row = ind(imin);
trial = LateGroup.TrialNumber(row);

% Fz examples: HHI8, trial 32 (AG) and trial 47 (AB)