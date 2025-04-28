
home_dir = '/project/3018051.01/ruggero/derivatives/first_level/task-rsa';
addpath('/project/3018051.01/ruggero/derivatives/');
cd(home_dir);

subjects = {'sub-102', 'sub-103', 'sub-104', 'sub-105', 'sub-106', 'sub-107', 'sub-108', 'sub-109', 'sub-110', 'sub-111', 'sub-112', 'sub-113', 'sub-114', 'sub-115', 'sub-116', 'sub-118', 'sub-119', 'sub-120', 'sub-121', 'sub-122', 'sub-123', 'sub-124', 'sub-125', 'sub-126', 'sub-127', 'sub-128', 'sub-129', 'sub-130', 'sub-131', 'sub-132', 'sub-133', 'sub-134', 'sub-135', 'sub-136', 'sub-201', 'sub-202', 'sub-203', 'sub-204', 'sub-205', 'sub-206', 'sub-207', 'sub-208', 'sub-209', 'sub-210', 'sub-211', 'sub-212', 'sub-213', 'sub-214', 'sub-215', 'sub-216', 'sub-218', 'sub-219', 'sub-220', 'sub-221', 'sub-222', 'sub-223', 'sub-224', 'sub-225', 'sub-226', 'sub-227', 'sub-228', 'sub-229', 'sub-230', 'sub-231', 'sub-232', 'sub-233', 'sub-234', 'sub-235'};
sessions = {'mri01', 'mri02'};

for k = 1:length(subjects)
    for j = 1:length(sessions)

        input_file = sprintf('/project/3018051.01/ruggero/derivatives/first_level/task-rsa/%s/%s_ses-%s_task-rsa_run-concat_events.tsv', subjects{k}, subjects{k}, sessions{j});
        output_file = sprintf('/project/3018051.01/ruggero/derivatives/first_level/task-rsa/%s/%s_SPM/%s_ses-%s_all-conditions.mat', subjects{k}, subjects{k}, subjects{k}, sessions{j});

        events = readtable(input_file, 'FileType', 'text', 'Delimiter', '\t');
        
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
        
        save(output_file, 'names', 'onsets', 'durations');

    end
end
