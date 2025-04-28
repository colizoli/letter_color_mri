
home_dir = '/project/3018051.01/ruggero/derivatives/first_level/task-rsa';
addpath('/project/3018051.01/ruggero/derivatives/');
cd(home_dir);

%eventpath = home_dir+"/session"+session+"/sub-"+this_subject+"/task-rsa_sub-"+this_subject+"_ses-0"+session+"_events.tsv";
events = readtable('/project/3018051.01/ruggero/derivatives/first_level/task-rsa/sub-202/sub-202_ses-mri02_task-rsa_run-concat_events.tsv', 'FileType', 'text', 'Delimiter', '\t');

% create oddbal condition
oddb_index = events.oddball == 1;
oddball_onsets = events(oddb_index,18);
duration_oddball = events(oddb_index,19);


% drops all rows in the events table that are oddball
g = events.oddball == 0;
events = events(g,:);

% onsets ALL GRAPHEMES (all letters except oddballs)
onsetGraphemes = events(:,18);

% get only COLORED GRAPHEMES
cg = strcmp(events.trial_type_color, 'color');
colored_graphemes = events(cg,:);
onsetColorGraphemes = colored_graphemes(:,18);
durationColorGraphemes = colored_graphemes(:,19);

% get only TRAINED GRAPHEMES
tg = strcmp(events.trial_type_letter, 'trained');
trained_graphemes = events(tg,:);
onsetTrainedGraphemes = trained_graphemes(:,18);
durationTrainedGraphemes = trained_graphemes(:,19);


%create all conditions
names = cell(1, 4);
onsets = cell(1, 4);
durations = cell(1, 4);

names{1} = 'all-graph';
names{2} = 'colored-graph';
names{3} = 'trained-graph';
names{4} = 'oddball';

onsets{1} = events(:,18);
onsets{2} = onsetColorGraphemes;
onsets{3} = onsetTrainedGraphemes;
onsets{4} = oddball_onsets;

durations{1} = events(:,19);
durations{2} = durationColorGraphemes;
durations{3} = durationTrainedGraphemes;
durations{4} = duration_oddball;

% 
for i = 1:length(onsets)
    onsets{i} = onsets{i}.onset;
    durations{i} = durations{i}.duration;
end


save('/project/3018051.01/ruggero/derivatives/first_level/task-rsa/sub-202/sub-202_SPM/sub-202_ses-mri02_all-conditions.mat', 'names', 'onsets', 'durations');