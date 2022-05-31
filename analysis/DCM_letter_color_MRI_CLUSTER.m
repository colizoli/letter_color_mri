

% linux cluster //
addpath('/project/2422045.01/faye/DCM/spm12/spm12/');
home_dir = '/project/2422045.01/faye/DCM';

cd(home_dir);
% spm fmri;

%subjects = {'201', '102', '202', '103', '203', '104','204','105','205','106','206','107','207','108','208','109','209','110','210','111','212','113','213','114','214','115','215','116','216','118','218','119','219','120','220','121','221','122','123','223','124','224','125','225','126', '226' };
subjects = {'201', '102', '202', '103', '203', '104','204','105','205','106','206','107','207','108','208'};


% 'normal' is 0000 up until and including 1985
nii_session1_1983 = {'223'};
nii_session1_1986 = {'213','214'};

% please check sub-106, is it 1985?
nii_session2_1983 = {'226'};
nii_session2_1986 = {'109','113','119','120','125','204'}


sessions = {'1'}; 

% load first level batch MAT TEMPLATE file
firstLevelSession = load("firstLevelSession_Template_CLUSTER.mat");
S1 = firstLevelSession.matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1);
% Note: to index structures, aim at indexing one level/field ABOVE the field you want to reach!
% S = 
% 
%   1×3 struct array with fields:
% 
%     scans
%     cond
%     multi
%     regress
%     multi_reg
%     hpf

% Copy subject 1's structure, then only cdphange what is necessary (onsets
% and file names)

for session_idx = 1:length(sessions)
    session = sessions(session_idx);
    
    for subject_idx = 1:length(subjects)        % subject loop
        this_subject = subjects{subject_idx};   % get current subject
        
        % path to events file for current subject and session
        this_events = home_dir+"/session"+session+"/sub-"+this_subject+"/task-rsa_sub-"+this_subject+"_ses-0"+session+"_events.tsv";
        this_nifti = home_dir+"/session"+session+"/sub-"+this_subject+"/task-rsa_sub-"+this_subject+"_ses-0"+session+"_bold_mni.tsv";
        % print out the file name
        disp(this_events);
        
        % calls the import data function and saves as a table called "events"
        events = import_events(this_events);
        
        % drops all rows in the events table that are oddball
        g = events.oddball == 0;
        events = events(g,:);
        
        % onsets ALL GRAPHEMES (all letters except oddballs)
        onsetGraphemes = events(:,16);
        
        % get only COLORED GRAPHEMES
        cg = events.trial_type_color == 'color';
        colored_graphemes = events(cg,:);
        onsetColorGraphemes = colored_graphemes(:,16);
        
        % get only TRAINED GRAPHEMES
        tg = events.trial_type_letter == 'trained';
        trained_graphemes = events(tg,:);
        onsetTrainedGraphemes = trained_graphemes(:,16);
        
        % get 1xNsubjects structure from firstLevelSession variable
        % copy the first subject's structure, then change only file names and onsets
        SS = S1;
        
        % graphemes
        SS.cond(1).onset = table2array(onsetGraphemes);
        % colored graphemes
        SS.cond(2).onset = table2array(onsetColorGraphemes);
        % trained graphemes
        SS.cond(3).onset = table2array(onsetTrainedGraphemes);
        
        % scans: {1986×1 cell}
        % loop through each volume file name, change subject number and session
        for nii = 1:length(SS.scans)
            %disp(nii);
            this_nii = SS.scans(nii);
            old = "sub-102";            % always first subject
            new = "sub-"+this_subject;  % current subject
            
            temp_string = strrep(this_nii{1},old,new);
            old = "session1";
            new = "session"+session;
            temp_string = strrep(temp_string,old,new);
            
            old = "ses-01";
            new = "ses-0"+session;
            
            SS.scans(nii) = cellstr(strrep(temp_string,old,new));
        end % end scan loop

        % change length of scans for some subjects
        if session_idx==1
            if any(strcmp(nii_session1_1983,this_subject));
                SS.scans = SS.scans(1:1984); % trim SS down
            elseif any(strcmp(nii_session1_1986,this_subject)); % add 1 extra scan
                temp_string = SS.scans(1986);
                old = '1985';
                new = '1986';
                SS.scans(1:1987)  = strrep(temp_string,old,new);
            end % end check scan length session 1
        elseif session_idx==2
            if any(strcmp(nii_session2_1983,this_subject));
                SS.scans = SS.scans(1:1984); % trim SS down
            elseif any(strcmp(nii_session2_1986,this_subject)); % add 1 extra scan
                temp_string = SS.scans(1986);
                old = '1985';
                new = '1986';
                SS.scans(1:1987)  = strrep(temp_string,old,new);
            end % end check scan length session 2
        end % end check session 

        % save entire structure for current subject
        firstLevelSession.matlabbatch{1, 1}.spm.stats.fmri_spec.sess(subject_idx) = SS;
    end % end subject loop
    
    % make sure to save the newly created firstLevelSession MAT file!!
    M_filename = "firstLevelSession"+session+".mat";
    matlabbatch = firstLevelSession.matlabbatch;
    save(M_filename, "matlabbatch");
    
end % end session loop


