#!/usr/bin/env python
# encoding: utf-8

"""
first_level.py
Created by O.Colizoli, June-2019
Last update: 14-01-2025
Python version 3.9

The following packages need to be installed because they are called from the command line, but not imported:
fsl

"""

# TO DO:

import os, subprocess, sys, glob
import shutil as sh
import nibabel as nib
import pandas as pd
import numpy as np
import json
from datetime import datetime
from IPython import embed as shell # for Oly's debugging only

class first_level_class(object):
    def __init__(self, subject, analysis_dir, bids_dir, deriv_dir, mask_dir, template_dir, timing_files_dir, TR):        
        self.subject        = 'sub-'+str(subject)
        self.analysis_dir   = str(analysis_dir)
        self.bids_dir       = str(bids_dir)
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.timing_files_dir = str(timing_files_dir)
        self.TR             = str(TR)
        
        # path to fmriprep output (preprocessed fmri data)
        self.fmriprep_dir = os.path.join(self.deriv_dir, 'fmriprep')
        
        # make directories for timing files, each task
        if not os.path.isdir(self.timing_files_dir):
            os.mkdir(self.timing_files_dir)  
        if not os.path.isdir(os.path.join(self.timing_files_dir, 'task-colors')):
            os.mkdir(os.path.join(self.timing_files_dir, 'task-colors'))
        if not os.path.isdir(os.path.join(self.timing_files_dir, 'task-letters')):
            os.mkdir(os.path.join(self.timing_files_dir, 'task-letters'))
        if not os.path.isdir(os.path.join(self.timing_files_dir, 'task-rsa')):
            os.mkdir(os.path.join(self.timing_files_dir, 'task-rsa'))
        
        # make directories for first level output (FSL), each task       
        self.first_level_dir = os.path.join(self.deriv_dir, 'first_level')
        if not os.path.isdir(self.first_level_dir):
            os.mkdir(self.first_level_dir)
        if not os.path.isdir(os.path.join(self.first_level_dir, 'task-colors')):
            os.mkdir(os.path.join(self.first_level_dir, 'task-colors'))
        if not os.path.isdir(os.path.join(self.first_level_dir, 'task-letters')):
            os.mkdir(os.path.join(self.first_level_dir, 'task-letters'))
        if not os.path.isdir(os.path.join(self.first_level_dir, 'task-rsa')):
            os.mkdir(os.path.join(self.first_level_dir, 'task-rsa'))
        
        # write unix commands to job to run in parallel
        job_dir = os.path.join(self.analysis_dir, 'jobs')
        if not os.path.isdir(job_dir):
            os.mkdir(job_dir)
        
        # write subject-specific job in analysis/jobs for bash scripts (overwrites on each call to class)
        self.first_level_job_path = os.path.join(self.analysis_dir, 'jobs', 'job_first_level_{}.txt'.format(self.subject))
        self.first_level_job = open(self.first_level_job_path, "w")
        self.first_level_job.write("#!/bin/bash\n")
        self.first_level_job.close()
            
    
    
    def loc_match_bold(self):
        """Match the LOC1 and LOC2 bold acquisition (nifti) files to the letters and colors localizer events. 
        """
        T = 6 # number of trailing timestamp characters in events file names
        
        localizer_df = pd.DataFrame()
        loc1 = []
        loc2 = []
        
        for session in ['ses-mri01', 'ses-mri02']:
            ###################
            # check acquisition time in JSON files
            ###################
            # JSON paths
            loc1_path = os.path.join(self.bids_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC01_dir-AP_bold.json'.format(self.subject, session))
            loc2_path = os.path.join(self.bids_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC02_dir-AP_bold.json'.format(self.subject, session))
        
            # LOC1 nifti
            with open(loc1_path) as f:
                loc1_json = json.load(f)
                loc1_acq = loc1_json['AcquisitionTime']
                loc1_acq = os.path.splitext(loc1_acq)[0] # drop extension
        
            # LOC2 nifti
            with open(loc2_path) as f:
                loc2_json = json.load(f)
                loc2_acq = loc2_json['AcquisitionTime']
                loc2_acq = os.path.splitext(loc2_acq)[0] # drop extension
            
            # letters events
            for f1 in os.listdir(os.path.join(self.bids_dir, self.subject, session, 'func')):
                if ('task-letters' in f1) and ('events' in f1):
                    fname, extension = os.path.splitext(f1) # split name from extension
                    letters_acq = fname[-T:]
        
            # colors events
            for f2 in os.listdir(os.path.join(self.bids_dir, self.subject, session, 'func')):
                if ('task-colors' in f2) and ('events' in f2):
                    fname, extension = os.path.splitext(f2) # split name from extension
                    colors_acq = fname[-T:]
        
            ### MATCHING ###
            letters_first = int(letters_acq) < int(colors_acq) # if letters are first, acquisition time is less
            loc1_first = datetime.strptime(loc1_acq, '%H:%M:%S') < datetime.strptime(loc2_acq, '%H:%M:%S') # if loc1 is first, acquisition time is less
            
            if letters_first:
                if loc1_first:
                    # letters == loc1
                    # colors = loc2
                    loc1.append('letters')
                    loc2.append('colors')                    
                else:
                    # colors = loc1
                    # letters  = loc2
                    loc1.append('colors')
                    loc2.append('letters')
        
            print('{} {} loc1_first={} letters_first={}'.format(self.subject, session, letters_first, loc1_first))       
        
        localizer_df['session'] = ['ses-mri01', 'ses-mri02']
        localizer_df['LOC01'] = loc1
        localizer_df['LOC02'] = loc2
        # save two copies in both localizer folders
        localizer_df.to_csv(os.path.join(self.first_level_dir, 'task-letters', '{}_localizer_matching.tsv'.format(self.subject)), sep='\t')
        localizer_df.to_csv(os.path.join(self.first_level_dir, 'task-colors', '{}_localizer_matching.tsv'.format(self.subject)), sep='\t')
        shell()
        print('success: loc_match_bold')
        
        
    def loc_combine_epi(self, task):
        """Concatenate the 2 sessions of EPI data to perform a single GLM on localizers.
        
        Args:
            task (str): which localizer task.
        
        Notes:
            The output is the concantenated bold of both sessions (input to first level).
            Also makes the session regressors as EV text files.
        """
        
        for preprocessed_tag in ['space-MNI152NLin6Asym_desc-preproc_bold', 'space-T1w_desc-preproc_bold']:
            
            df_trs = pd.DataFrame() # save number of TRs for events file concatenation
        
            # output is the concantenated bold of all runs per session (input to first level)
            out_file = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_ses-concat_task-{}_{}.nii.gz'.format(self.subject, task, preprocessed_tag))
        
            for s, session in enumerate(['ses-mri01', 'ses-mri02']): 
                for this_file in os.listdir(os.path.join(self.fmriprep_dir, self.subject, session, 'func')):
                    shell()
                    if (task in this_file) and (session in this_file) and (preprocessed_tag in this_file) and ('nii.gz' in this_file):
                        
                        shell()
                        nii = nib.load(nii_path)
                        # count TRs
                        df_trs[this_run] = [nii.get_fdata().shape[-1]] # 4th dimension
                        # append to list
                        bold.append(nii.get_fdata())
                        if '1' in this_run:
                            save_path = nii_path
                        
            bold = np.concatenate(bold, axis=-1) # concatenate list to one long ndarray
        
            # get nifti header from first run
            n1 = nib.load(save_path)
            out_data = nib.Nifti1Image(bold, affine=n1.affine, header=n1.header) # pass affine and header from last MNI image
            out_data.set_data_dtype(np.float32)
            nib.save(out_data, out_file)
        
            # save trs per session
            df_trs.to_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_TRs.tsv'.format(self.subject, session, task)), sep='\t')
            print(bold.shape)
            
        
            # output session regressors as custom 1 column EV files
            ses1 = np.concatenate( (np.repeat(1,len(bold1)),  np.repeat(0,len(bold2)) ), axis=0) # session 1: 1s run 1
            ses2 = np.concatenate( (np.repeat(0,len(bold1)),  np.repeat(1,len(bold2)) ), axis=0) # session 2: 1s run 2
            # ses-01 EV
            out_ses = os.path.join(self.deriv_dir, 'timing_files','task-{}'.format(task), 'task-{}_{}_{}.txt'.format(task, self.subject, 'ses-01'))
            output = np.array(ses1.T) # 1 x TR
            np.savetxt(out_fn, output, delimiter='/t', fmt='%i') #column has to be an integer for FSL!
            # ses-02 EV
            out_ses = os.path.join(self.deriv_dir, 'timing_files','task-{}'.format(task), 'task-{}_{}_{}.txt'.format(task, self.subject, 'ses-02'))
            output = np.array(ses2.T) # 1 x TR
            np.savetxt(out_fn, output, delimiter='/t', fmt='%i') #column has to be an integer for FSL!

            print('success: loc_combine_epi {}'.format(self.subject))
     
     
    def loc_combine_timing_files(self, task):
        """Concatenate the timing files of 2 sessions to perform a single GLM on localizers.
        
        Args:
            task (str): which localizer task.
        
        Notes:
            2nd session have to add time = #TRs * TR
            GLM timing files for ABAB blocked design (localizers)
        """
        STIM_DUR = 0.75 # stimulus duration seconds
        
        n1_fn = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, 'ses-01', task), 'filtered_func_data.nii.gz')
        n2_fn = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, 'ses-02', task), 'filtered_func_data.nii.gz')
        
        # check if both localizers present
        if os.path.exists(n1_fn) and os.path.exists(n2_fn):
            
            # take FIRST session's BOLD to count TRs to add to 2nd sessions' onsets
            bold1 = nib.load(n1_fn)
            ntrs = bold1.shape[-1]  # number of TRs
            time2add = ntrs*float(self.TR) # time to add in seconds to block 2's onsets
            print('Time to add to run 2: {}'.format(time2add))
        
            # open block 2's events
            events2 = pd.read_csv(os.path.join(self.deriv_dir, self.subject, 'ses-02', 'func', '{}_{}_task-{}_events.tsv'.format(self.subject, 'ses-02', task)), sep='\t')
            events2['onset'] =  events2['onset'] + time2add
            
            # open block 1's events and concantenate
            events1 = pd.read_csv(os.path.join(self.deriv_dir, self.subject, 'ses-01', 'func', '{}_{}_task-{}_events.tsv'.format(self.subject, 'ses-01', task)), sep='\t')
            events = pd.concat([events1, events2], axis=0)
            events.to_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), 'task-{}_{}_events.tsv'.format(task, self.subject)), sep='\t')  # save concantenated events file

            # generate 3 column files for each of the 2x2 conditions
            for c,cond in enumerate(np.unique(events['trial_type'])):
                out_fn = os.path.join(self.deriv_dir, 'timing_files','task-{}'.format(task), 'task-{}_{}_{}.txt'.format(task,self.subject,cond))

                # main regressors
                first = np.array(events[events['trial_type']==cond]['onset']) # onset in s
                second = np.repeat(STIM_DUR, len(first))    # duration in s
                third = np.array(np.repeat(1, len(first)),dtype=int)    # amplitude
                output = np.array(np.vstack((first, second, third)).T) # 1 x 3
                np.savetxt(out_fn, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
                print(out_fn)
        else:
            print('{} {} missing task-{}!'.format(self.subject, session, task))
        print('success: loc_combine_timing_files {}'.format(self.subject))
    
    
    def loc_combine_motion_regressors(self, task):
        """Concatenate the 2 sessions of the motion regressors to perform a single GLM.
        
        Args:
            task (str): which localizer task.
        
        Notes:
            The output is the concantenated motion regressors as nifti files of both sessions.
        """
        
        #### Motion parameters ####
        n1_fn = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, 'ses-01', task), 'filtered_func_data.nii.gz')
        n2_fn = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, 'ses-02', task), 'filtered_func_data.nii.gz')
        
        # check if both localizers present
        if os.path.exists(n1_fn) and os.path.exists(n2_fn):
            # 6 motion regressors
            for mcf in np.arange(6):
                mcf_out = os.path.join(self.first_level_dir, 'task-{}'.format(task), 'task-{}_{}_mcf_par{}.nii.gz'.format(task, self.subject. mcf + 1)) 
                
                mcf_fn1 = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, 'ses-01', task), 'mc', 'mcf_par{}.nii.gz'.format(mcf + 1))
                mcf_fn2 = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, 'ses-02', task), 'mc', 'mcf_par{}.nii.gz'.format(mcf + 1))
                n1 = nib.load(mcf_fn1) # preprocessed session 1
                n2 = nib.load(mcf_fn2) # preprocessed session 2
            
                bold1 = n1.get_data()
                bold2 = n2.get_data()
                bold = np.concatenate([bold1, bold2],axis=-1)
            
                nii_out = nib.Nifti1Image(bold, affine=n1.affine, header=n1.header) # pass affine and header from last MNI image
                nii_out.set_data_dtype(np.float32)
                nib.save(nii_out, mri_out)        
        else:
            print('{} {} missing task-{}!'.format(self.subject, session, task))
           
        print('success: loc_combine_motion_regressors {}'.format(self.subject))
        

    def loc_fsf(self,task):
        # Creates the FSF files for each subject's first level analysis - localizers
        # Run the actual FSF from the command line: feat task-colors_sub-01_ses-01.fsf
            
        template_filename = os.path.join(self.analysis_dir,'templates','task-{}_first_level_template.fsf'.format(task))
    
        markers = [
            '[$OUTPUT_PATH]', 
            '[$NR_TRS]', 
            '[$INPUT_FILENAME]', 
            '[$NUISANCE]', 
            '[$EV1_FILENAME]',
            '[$EV2_FILENAME]', 
            '[$NR_VOXELS]',
            '[$MNI_BRAIN]'
        ]
        
        FSF_filename = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}.fsf'.format(task,self.subject)) # save fsf
        output_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}'.format(task,self.subject)) 
        
        BOLD = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_bold_mni.nii.gz'.format(task,self.subject)) 
        # calculate size of input data
        nii = nib.load(BOLD).get_data() # only do once 
        nr_trs = str(nii.shape[-1])
        nr_voxels = str(nii.size)
        
        # motion parameters and columns for each session's mean
        nuisance_regressors = os.path.join(self.timing_files_dir,'task-{}'.format(task),'task-{}_{}_nuisance_regressors.txt'.format(task,self.subject))
        
        # timing files for each EV
        if task == 'colors':
            EVS = ['Color','Black']
        elif task == 'letters':
            EVS = ['Letter','Symbol']
        EV1_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}.txt'.format(task,self.subject,EVS[0]))
        EV2_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}.txt'.format(task,self.subject,EVS[1]))
        
        MNI_BRAIN  = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz')
        
        # replacements
        replacements = [ # needs to match order of 'markers'
            output_path,
            nr_trs,
            BOLD,
            nuisance_regressors,
            EV1_path,
            EV2_path,
            nr_voxels,
            MNI_BRAIN
        ]

        # open the template file, load the text data
        f = open(template_filename,'r')
        filedata = f.read()
        f.close()

        # search and replace
        for st,this_string in enumerate(markers):
            filedata = filedata.replace(this_string,replacements[st])

        # write output file
        f = open(FSF_filename,'w')
        f.write(filedata)
        f.close()
    
        # open job and write command as new line
        cmd = 'feat {}'.format(FSF_filename)
        self.first_level_job = open(self.first_level_job_path, "a") # append is important, not write
        self.first_level_job.write(cmd)   # feat command
        self.first_level_job.write("\n\n")  # new line
        self.first_level_job.close()
        print('success: loc_fsf {}'.format(FSF_filename))
        
    
    def rsa_combine_epi(self, task='rsa'):
        """Concatenate the 4 runs per session of EPI data to perform a single GLM (RSA task).
        
        Args:
            task (str): which task. Default 'rsa'
        
        Notes:
            The output is the concantenated bold of all runs (input to first level).
            TO DO: Equalize number of runs per session per participant (in case of missing runs).
        """
                
        for preprocessed_tag in ['space-MNI152NLin6Asym_desc-preproc_bold', 'space-T1w_desc-preproc_bold']:
        
            for session in ['ses-mri01','ses-mri02']:
            
                df_trs = pd.DataFrame() # save number of TRs for events file concatenation
            
                # output is the concantenated bold of all runs per session (input to first level)
                out_file = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_{}.nii.gz'.format(self.subject, session, task, preprocessed_tag))
            
                bold = []
                for r, this_run in enumerate(['run-1', 'run-2', 'run-3', 'run-4']): 
                    for this_file in os.listdir(os.path.join(self.fmriprep_dir, self.subject, session, 'func')):
                        if (this_run in this_file) and (preprocessed_tag in this_file) and ('nii.gz' in this_file):
                            nii_path = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500RSA_dir-AP_{}_{}.nii.gz'.format(self.subject, session, this_run, preprocessed_tag))
                            nii = nib.load(nii_path)
                            # count TRs
                            df_trs[this_run] = [nii.get_fdata().shape[-1]] # 4th dimension
                            # append to list
                            bold.append(nii.get_fdata())
                            if '1' in this_run:
                                save_path = nii_path
                            
                bold = np.concatenate(bold, axis=-1) # concatenate list to one long ndarray
            
                # get nifti header from first run
                n1 = nib.load(save_path)
                out_data = nib.Nifti1Image(bold, affine=n1.affine, header=n1.header) # pass affine and header from last MNI image
                out_data.set_data_dtype(np.float32)
                nib.save(out_data, out_file)
            
                # save trs per session
                df_trs.to_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_TRs.tsv'.format(self.subject, session, task)), sep='\t')
                print(bold.shape)
                    
        print('success: rsa_combine_epi')
        
        
    def rsa_combine_events(self, task='rsa'):
        """For the RSA task, concantenate the events files of all runs and output in first_level directory.
        
        Args:
            task (str): which task. Default 'rsa'
        
        Notes:
            To DO Equalize number of runs per session per participant (in case of missing runs).
        """
        
        preprocessed_tag = 'space-T1w_desc-preproc_bold'
        
        for session in ['ses-mri01','ses-mri02']:
            
            # open tr count from rsa_combine_epi()
            df_trs = pd.read_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_TRs.tsv'.format(self.subject, session, task)), sep='\t')
            df_trs = df_trs.loc[:, ~df_trs.columns.str.contains('^Unnamed')]
                        
            rsa_runs = df_trs.columns.values
            nruns = len(df_trs.columns.values)
            
            time2add = [0] # don't add anytime to run 1
            
            # check how long each run was and calculate seconds
            for this_run in rsa_runs:
                ntrs = df_trs[this_run][0]  # number of TRs
                time2add.append(ntrs*float(self.TR)) # time to add in seconds to next run's onsets
                print('Time to add to {}: {}'.format(this_run, ntrs*float(self.TR)))
                
            time2add = np.cumsum(time2add) # cummulative sum per run
            
            # open events files and add times
            all_events = pd.DataFrame()
            for r, this_run in enumerate(rsa_runs):
                for this_file in os.listdir(os.path.join(self.bids_dir, self.subject, session, 'func')):
                    if (task in this_file) and (this_run in this_file) and ('events' in this_file):
                        print(this_file)
                        print(this_run)
                        print(time2add[r])
                        events = pd.read_csv(os.path.join(self.bids_dir, self.subject, session, 'func', this_file), sep='\t')
                        events['onset'] =  events['onset'] + time2add[r]
                        all_events = pd.concat([all_events, events], axis=0)
                print(all_events.shape)
            
            # shell()
            # add unique identifiers for each color
            rgb_codes = [
                (all_events['r'] == 188) & (all_events['g'] == 188) & (all_events['b'] == 188), # grey (oddballs)
                (all_events['r'] == 117) & (all_events['g'] == 117) & (all_events['b'] == 117), # grey (oddballs)
                (all_events['r'] == 128) & (all_events['g'] == 128) & (all_events['b'] == 128), # grey (oddballs)
                (all_events['r'] == 0) & (all_events['g'] == 0) & (all_events['b'] == 0), # black
                (all_events['r'] == 0) & (all_events['g'] == 163) & (all_events['b'] == 228), # light_blue
                (all_events['r'] == 161) & (all_events['g'] == 199) & (all_events['b'] == 70), # lime_green
                (all_events['r'] == 183) & (all_events['g'] == 61) & (all_events['b'] == 160), # magenta
                (all_events['r'] == 181) & (all_events['g'] == 44) & (all_events['b'] == 67), # dark_red
                (all_events['r'] == 16) & (all_events['g'] == 114) & (all_events['b'] == 86), # dark_green
                (all_events['r'] == 237) & (all_events['g'] == 114) & (all_events['b'] == 162), # pink
                (all_events['r'] == 58) & (all_events['g'] == 175) & (all_events['b'] == 75), # green
                (all_events['r'] == 248) & (all_events['g'] == 154) & (all_events['b'] == 28), # light_orange
                (all_events['r'] == 109) & (all_events['g'] == 57) & (all_events['b'] == 142), # purple
                (all_events['r'] == 239) & (all_events['g'] == 79) & (all_events['b'] == 41), # orange
                (all_events['r'] == 49) & (all_events['g'] == 60) & (all_events['b'] == 163), # blue
                (all_events['r'] == 255) & (all_events['g'] == 211) & (all_events['b'] == 0), # yellow
                (all_events['r'] == 9) & (all_events['g'] == 181) & (all_events['b'] == 172) # teal
            ]
            
            color_names = [
                'grey',
                'grey',
                'grey',
                'black',
                'light_blue',
                'lime_green',
                'magenta',
                'dark_red',
                'dark_green',
                'pink',
                'green',
                'light_orange',
                'purple',
                'orange',
                'blue',
                'yellow',
                'teal'
            ]
            all_events['color_name'] = np.select(rgb_codes, color_names)
                    
            # save concantenated events file
            all_events.to_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_events.tsv'.format(self.subject, session, task)), sep='\t') 

        print('success: rsa_combine_events')    


    def rsa_dcm_split_nifti(self, task='rsa'):
        """Split the concatenated four runs of the RSA task into single images for SPM (DCM analysis) and unzip compressed niftis.
        
        Args:
            task (str): which task. Default 'rsa'.
        
        Notes:
            See rsa_combine_epi()
        """
        
        preprocessed_tag = 'space-T1w_desc-preproc_bold'
        
        for session in ['ses-mri01','ses-mri02']:            
        
            # input is the concantenated bold of all runs per session (input to first level)
            nii_in = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_{}.nii.gz'.format(self.subject, session, task, preprocessed_tag))
            
            nii_out_base = os.path.join(self.deriv_dir, 'DCM', self.subject, session, '{}_{}_task-{}_run-concat_'.format(self.subject, session, task))
            
            # make subject/session folders inside derivatives/DCM if it doesn't exist
            path_subj = os.path.join(self.deriv_dir, 'DCM', self.subject)
            if not os.path.isdir(path_subj):
                os.mkdir(path_subj)
            path_subj_sess = os.path.join(path_subj, session)
            if not os.path.isdir(path_subj_sess):
                os.mkdir(path_subj_sess)    
                                
            # fslsplit input output_base_name -t
            cmd = 'fslsplit {} {} -t'.format(nii_in, nii_out_base)
            print(cmd)
            results = subprocess.call(cmd, shell=True, bufsize=0)
            
            # unzip the compressed nifti files for SPM, overwrites original .gz files
            for nii in os.listdir(path_subj_sess):
                file_path = os.path.join(path_subj_sess, nii)
                cmd = 'gzip -d {}'.format(file_path)
                print(cmd)
                results = subprocess.call(cmd, shell=True, bufsize=0)

        print('success: rsa_dcm_split_nifti')    
        
        
    def rsa_nuisance_regressors(self, task='rsa'):
        """Combine all TR-based nuisance regressors for all runs into a single file.
        
        Args:
            task (str): which task. Default 'rsa'.
        
        Notes:
            Includes the following nuisance regressors:
            1. first 5 principle components of the compcor output of fmriprep (white matter and csf).
            2. 6 motion parameters from fmriprep.
            3. cosine regressors for low-frequency drift.
            4. the 4 run regressors (for concantenations).
                    
            confounds file name: derivatives/fmriprep/sub-xxx/ses-mri-xx/func/sub-xxx_ses-mrixx_task-cmrr2isomb4TR1500RSA_dir-AP_run-x_desc-confounds_timeseries.tsv
            anatomical compcor column names in confounds file: a_comp_cor_00, a_comp_cor_01, a_comp_cor_02, a_comp_cor_03, a_comp_cor_04
            motion column names in confounds file: trans_x, trans_y, trans_z, rot_x, rot_y, rot_z
        """
        
        for session in ['ses-mri01','ses-mri02']:
            
            save_regressors = pd.DataFrame() # output file
            
            # open tr count from rsa_combine_epi()
            df_trs = pd.read_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_TRs.tsv'.format(self.subject, session, task)), sep='\t')
            df_trs = df_trs.loc[:, ~df_trs.columns.str.contains('^Unnamed')]
                        
            rsa_runs = df_trs.columns.values # use these runs
            nruns = len(df_trs.columns.values) # number of trs per run
            
            for r in rsa_runs:
                this_run_regressors = pd.DataFrame()
                # output of fmriprep per subject, per session, per run
                confounds = pd.read_csv(os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500RSA_dir-AP_{}_desc-confounds_timeseries.tsv'.format(self.subject, session, r)), sep='\t', float_precision='high')
                # compcor
                this_run_regressors['a_comp_cor_00'] = np.array(confounds['a_comp_cor_00'])
                this_run_regressors['a_comp_cor_01'] = np.array(confounds['a_comp_cor_01'])
                this_run_regressors['a_comp_cor_02'] = np.array(confounds['a_comp_cor_02'])
                this_run_regressors['a_comp_cor_03'] = np.array(confounds['a_comp_cor_03'])
                this_run_regressors['a_comp_cor_04'] = np.array(confounds['a_comp_cor_04'])
                # motion
                this_run_regressors['trans_x'] = np.array(confounds['trans_x'])
                this_run_regressors['trans_y'] = np.array(confounds['trans_y'])
                this_run_regressors['trans_z'] = np.array(confounds['trans_z'])
                this_run_regressors['rot_x'] = np.array(confounds['rot_x'])
                this_run_regressors['rot_y'] = np.array(confounds['rot_y'])
                this_run_regressors['rot_z'] = np.array(confounds['rot_z'])
                # cosine (low-frequency drift)
                this_run_regressors['cosine00'] = np.array(confounds['cosine00'])
                this_run_regressors['cosine01'] = np.array(confounds['cosine01'])
                this_run_regressors['cosine02'] = np.array(confounds['cosine02'])
                this_run_regressors['cosine03'] = np.array(confounds['cosine03'])
                this_run_regressors['cosine04'] = np.array(confounds['cosine04'])
                this_run_regressors['cosine05'] = np.array(confounds['cosine05'])
                this_run_regressors['cosine06'] = np.array(confounds['cosine06'])
                this_run_regressors['cosine07'] = np.array(confounds['cosine07'])
                this_run_regressors['cosine08'] = np.array(confounds['cosine08'])
                this_run_regressors['cosine09'] = np.array(confounds['cosine09'])
                # runs
                this_run_regressors[r] = np.repeat(1, confounds.shape[0])
                
                save_regressors = pd.concat([save_regressors, this_run_regressors], axis=0)
            
            save_regressors = save_regressors.fillna(0)   
            # save without headers and without index, tab separated
            save_regressors.to_csv(
                os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_nuisance_regressors.txt'.format(self.subject, session, task)),
                sep='\t', float_format='%.16f', header=None, index_label=False, index=False)
            
        print('success: rsa_nuisance_regressors {}'.format(self.subject))
    

    def rsa_timing_files_oddballs(self, task='rsa'):
        """Create GLM timing files for RSA - oddball trials.
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            Event-related design.
            Duration is stimulus duration, which is preferred because stimulus was usually longer than 
            RT and sometimes people missed button presses or technical errors with recording.
        """

        for session in ['ses-mri01','ses-mri02']: # load session data with onsets
            events = pd.read_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_events.tsv'.format(self.subject, session, task)), sep='\t') # save concantenated events file
            
            #######################
            ### odd ball trials ###
            out_file = os.path.join(self.timing_files_dir, 'task-rsa', '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'oddballs'))
            # main regressors
            first = np.array(events[(events['oddball']==1)]['onset']) # onset in s
            second = events[(events['oddball']==1)]['duration'] # duration in s
            third = np.array(np.repeat(1, len(first)), dtype=int) # amplitude
            output = np.array(np.vstack((first, np.array(second), third)).T) # 1 x 3
            np.savetxt(out_file, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
            print(out_file)
            
        print('success: rsa_timing_files_oddballs')    
        

    def rsa_timing_files_letters(self, task='rsa'):
        """Create GLM timing files for RSA - each letter in each color (2 trials per run).
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            Event-related design.
            Outputs timing file for the oddball trials as well (nuisance EV).
        """

        color_name = [
                # 'grey',
                # 'grey',
                # 'grey',
                'black',
                'light_blue',
                'lime_green',
                'magenta',
                'dark_red',
                'dark_green',
                'pink',
                'green',
                'light_orange',
                'purple',
                'orange',
                'blue',
                'yellow',
                'teal'
            ]

        for session in ['ses-mri01','ses-mri02']: # load session data with onsets
            events = pd.read_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_events.tsv'.format(self.subject, session, task)), sep='\t') # save concantenated events file


            # generate 3 column files for each of the conditions
            for l,lcond in enumerate(np.unique(events['letter'])): # letters
                for c,ccond in enumerate(color_name): # 13 trained colors
                    # this letter in this color?
                    this_ev = events[(events['letter']==lcond) & (events['color_name']==ccond)]
                    # check if letter-color condition exists to prevent writing out empty EV timing files
                    if not this_ev.empty: 
                        # all files have to have the same names for subjects
                        out_file = os.path.join(self.timing_files_dir, 'task-rsa', '{}_{}_task-{}_{}_{}.txt'.format(self.subject, session, task, lcond, np.array(this_ev['trial_type_color'])[0]))
                        # main regressors
                        first = np.array(this_ev['onset']) # onset in s
                        second = np.array(this_ev['duration']) # onset in s # duration in s
                        third = np.array(np.repeat(1, len(first)),dtype=int) # amplitude
                        output = np.array(np.vstack((first, second, third)).T) # 1 x 3
                        np.savetxt(out_file, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
                        print(out_file)
        print('success: rsa_timing_files_letters')    
        
    
    def rsa_timing_files_2x2(self, task='rsa'):
        """Create GLM timing files for RSA - simple 2x2 comparison: trained/untrained vs. color/black
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            Event-related design.
        """
        
        for session in ['ses-mri01','ses-mri02']:
            events = pd.read_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_events.tsv'.format(self.subject, session, task)), sep='\t') # save concantenated events file
            
            # generate 3 column files for each of the 2x2 conditions
            for l,lcond in enumerate(['trained', 'untrained']): # letter condition
                for c,ccond in enumerate(['black', 'color']): # color condition
        
                    out_file = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_{}.txt'.format(self.subject, session, task, lcond, ccond))
                    # main regressors
                    first = np.array(events[(events['trial_type_letter']==lcond) & (events['trial_type_color']==ccond)]['onset']) # onset in s
                    second = np.array(events[(events['trial_type_letter']==lcond) & (events['trial_type_color']==ccond)]['duration']) # duration in s
                    third = np.array(np.repeat(1, len(first)),dtype=int) # amplitude
                    output = np.array(np.vstack((first, second, third)).T) # 1 x 3
                    np.savetxt(out_file, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
                    print(out_file)
        print('success: rsa_timing_files_2x2')    
        
    
    def rsa_letters_fsf(self, task='rsa'):
        """Create the FSF files for each subject's first level analysis - RSA design each letter in each color.
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            Subject native space.
            Event-related design.
            Oddball trials are nuisance EV, physiological and motion regressors are listed in the confound text file.
            Run the actual FSF as batch script or from the command line: feat task-rsa_sub-01_ses-01.fsf
        """
        
        preprocessed_tag = 'space-T1w_desc-preproc_bold'
        
        template_filename = os.path.join(self.template_dir, 'task-{}_letters_first_level_template.fsf'.format(task))

        markers = [
            '[$OUTPUT_PATH]',
            '[$NR_TRS]',
            '[$NR_VOXELS]',
            '[$INPUT_FILENAME]',
            '[$NUISANCE]',
            '[$EV1_FILENAME]','[$EV2_FILENAME]','[$EV3_FILENAME]','[$EV4_FILENAME]','[$EV5_FILENAME]','[$EV6_FILENAME]','[$EV7_FILENAME]','[$EV8_FILENAME]','[$EV9_FILENAME]', '[$EV10_FILENAME]',
            '[$EV11_FILENAME]','[$EV12_FILENAME]','[$EV13_FILENAME]','[$EV14_FILENAME]','[$EV15_FILENAME]','[$EV16_FILENAME]','[$EV17_FILENAME]','[$EV18_FILENAME]','[$EV19_FILENAME]', '[$EV20_FILENAME]',
            '[$EV21_FILENAME]','[$EV22_FILENAME]','[$EV23_FILENAME]','[$EV24_FILENAME]','[$EV25_FILENAME]','[$EV26_FILENAME]','[$EV27_FILENAME]','[$EV28_FILENAME]','[$EV29_FILENAME]', '[$EV30_FILENAME]',
            '[$EV31_FILENAME]','[$EV32_FILENAME]','[$EV33_FILENAME]','[$EV34_FILENAME]','[$EV35_FILENAME]','[$EV36_FILENAME]','[$EV37_FILENAME]','[$EV38_FILENAME]','[$EV39_FILENAME]', '[$EV40_FILENAME]',
            '[$EV41_FILENAME]','[$EV42_FILENAME]','[$EV43_FILENAME]','[$EV44_FILENAME]','[$EV45_FILENAME]','[$EV46_FILENAME]','[$EV47_FILENAME]','[$EV48_FILENAME]','[$EV49_FILENAME]', '[$EV50_FILENAME]',
            '[$EV51_FILENAME]','[$EV52_FILENAME]',
            '[$EV53_FILENAME]' # oddballs
        ]

        for session in ['ses-mri01','ses-mri02']:
            fsf_filename = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_letters.fsf'.format(self.subject, session, task)) # save fsf
            output_path = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_letters'.format(self.subject, session, task))

            bold = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_{}.nii.gz'.format(self.subject, session, task, preprocessed_tag)) 
            # calculate size of input data
            nii = nib.load(bold).get_fdata() # only do once
            nr_trs = str(nii.shape[-1])
            nr_voxels = str(nii.size)

            # nuiscance regressors
            nuisance_regressors = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_nuisance_regressors.txt'.format(self.subject, session, task))

            # timing files for each EV
            # black
            EV1_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'a_black'))
            EV2_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'b_black'))
            EV3_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'c_black'))
            EV4_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'d_black'))
            EV5_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'e_black'))
            EV6_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'f_black'))
            EV7_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'g_black'))
            EV8_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'h_black'))
            EV9_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'i_black'))
            EV10_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'j_black'))
            EV11_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'k_black'))
            EV12_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'l_black'))
            EV13_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'m_black'))
            EV14_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'n_black'))
            EV15_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'o_black'))
            EV16_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'p_black'))
            EV17_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'q_black'))
            EV18_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'r_black'))
            EV19_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'s_black'))
            EV20_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'t_black'))
            EV21_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'u_black'))
            EV22_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'v_black'))
            EV23_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'w_black'))
            EV24_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'x_black'))
            EV25_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'y_black'))
            EV26_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'z_black'))
            # color
            EV27_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'a_color'))
            EV28_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'b_color'))
            EV29_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'c_color'))
            EV30_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'d_color'))
            EV31_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'e_color'))
            EV32_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'f_color'))
            EV33_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'g_color'))
            EV34_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'h_color'))
            EV35_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'i_color'))
            EV36_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'j_color'))
            EV37_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'k_color'))
            EV38_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'l_color'))
            EV39_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'m_color'))
            EV40_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'n_color'))
            EV41_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'o_color'))
            EV42_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'p_color'))
            EV43_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'q_color'))
            EV44_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'r_color'))
            EV45_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'s_color'))
            EV46_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'t_color'))
            EV47_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'u_color'))
            EV48_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'v_color'))
            EV49_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'w_color'))
            EV50_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'x_color'))
            EV51_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'y_color'))
            EV52_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task,'z_color'))
            # oddballs last
            EV53_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'oddballs'))

            # replacements
            replacements = [ # needs to match order of 'markers'
                output_path,
                nr_trs,
                nr_voxels,
                bold,
                nuisance_regressors,
                EV1_path, EV2_path, EV3_path, EV4_path, EV5_path, EV6_path, EV7_path, EV8_path, EV9_path, EV10_path,
                EV11_path, EV12_path, EV13_path, EV14_path, EV15_path, EV16_path, EV17_path, EV18_path, EV19_path, EV20_path,
                EV21_path, EV22_path, EV23_path, EV24_path, EV25_path, EV26_path, EV27_path, EV28_path, EV29_path, EV30_path,
                EV31_path, EV32_path, EV33_path, EV34_path, EV35_path, EV36_path, EV37_path, EV38_path, EV39_path, EV40_path,
                EV41_path, EV42_path, EV43_path, EV44_path, EV45_path, EV46_path, EV47_path, EV48_path, EV49_path, EV50_path,
                EV51_path, EV52_path,
                EV53_path
            ]

            # open the template file, load the text data
            f = open(template_filename,'r')
            filedata = f.read()
            f.close()

            # search and replace
            for st,this_string in enumerate(markers):
                filedata = filedata.replace(this_string, replacements[st])

            # write output file
            f = open(fsf_filename,'w')
            f.write(filedata)
            f.close()
    
            # open preprocessing job and write command 
            cmd = 'feat {}'.format(fsf_filename)
            self.first_level_job = open(self.first_level_job_path, "a") # append is important, not write
            self.first_level_job.write(cmd)   # feat command
            self.first_level_job.write("\n\n")  # new line
            self.first_level_job.close()
        print('success: rsa_letters_fsf {}'.format(fsf_filename))
    
    
    def rsa_2x2_fsf(self, task='rsa'):
        """Createf the FSF files for each subject's first level analysis - RSA 2x2 trained/untrained vs. black/color
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            MNI Space
            Event-related design.
            Oddball trials are nuisance EV, physiological and motion regressors are listed in the confound text file.
            Run the actual FSF as batch script or from the command line: feat task-rsa_sub-01_ses-01.fsf
        """
        preprocessed_tag = 'space-MNI152NLin6Asym_desc-preproc_bold'
        
        template_filename = os.path.join(self.template_dir, 'task-rsa_2x2_first_level_template.fsf')
    
        markers = [
            '[$OUTPUT_PATH]', 
            '[$NR_TRS]', 
            '[$INPUT_FILENAME]', 
            '[$NUISANCE]', 
            '[$EV1_FILENAME]',
            '[$EV2_FILENAME]', 
            '[$EV3_FILENAME]', 
            '[$EV4_FILENAME]', 
            '[$EV5_FILENAME]', 
            '[$NR_VOXELS]',
        ]
        
        for session in ['ses-mri01','ses-mri02']:
            fsf_filename = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_2x2.fsf'.format(self.subject, session, task)) # save fsf
            output_path = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_2x2'.format(self.subject, session, task)) 
        
            bold = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_{}.nii.gz'.format(self.subject, session, task, preprocessed_tag)) 
            # calculate size of input data
            nii = nib.load(bold).get_fdata() # only do once 
            nr_trs = str(nii.shape[-1])
            nr_voxels = str(nii.size)
        
            # nuisance regressors
            nuisance_regressors = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_nuisance_regressors.txt'.format(self.subject, session, task))
        
            # timing files for each EV
            EV1_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'trained_black'))
            EV2_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'trained_color'))
            EV3_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'untrained_black'))
            EV4_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'untrained_color'))
            EV5_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, session, task, 'oddballs'))
                
            # replacements
            replacements = [ # needs to match order of 'markers'
                output_path,
                nr_trs,
                bold,
                nuisance_regressors,
                EV1_path,
                EV2_path,
                EV3_path,
                EV4_path,
                EV5_path,
                nr_voxels,
            ]

            # open the template file, load the text data
            f = open(template_filename,'r')
            filedata = f.read()
            f.close()

            # search and replace
            for st,this_string in enumerate(markers):
                filedata = filedata.replace(this_string, replacements[st])

            # write output file
            f = open(fsf_filename,'w')
            f.write(filedata)
            f.close()
    
            # open preprocessing job and write command as new line
            cmd = 'feat {}'.format(fsf_filename)
            self.first_level_job = open(self.first_level_job_path, "a") # append is important, not write
            self.first_level_job.write(cmd)   # feat command
            self.first_level_job.write("\n\n")  # new line
            self.first_level_job.close()
        print('success: rsa_2x2_fsf {}'.format(fsf_filename))


