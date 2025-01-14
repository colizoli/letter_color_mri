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
            os.mkdir(os.path.join(self.timing_files_dir, 'task-colors'))
            os.mkdir(os.path.join(self.timing_files_dir, 'task-letters'))
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
            
        self.first_level_job_path = os.path.join(self.analysis_dir, 'jobs', 'job_first_level_{}.txt'.format(self.subject))
        if not os.path.exists(self.first_level_job_path):
            self.first_level_job = open(self.first_level_job_path, "w")
            self.first_level_job.write("#!/bin/bash\n")
            self.first_level_job.close()
            
            
    def loc_combine_epi(self, task):
        """Concatenate the 2 sessions of EPI data to perform a single GLM on localizers.
        
        Args:
            task (str): which localizer task.
        
        Notes:
            The output is the concantenated bold of both sessions (input to first level).
            Also makes the session regressors as EV text files.
        """
        mri_out = os.path.join(self.first_level_dir, 'task-{}'.format(task), 'task-{}_{}_bold.nii.gz'.format(task, self.subject)) 
                
        n1_fn = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, 'ses-01', task), 'filtered_func_data.nii.gz')
        n2_fn = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, 'ses-02', task), 'filtered_func_data.nii.gz')
        
        # check if both localizers present
        if os.path.exists(n1_fn) and os.path.exists(n2_fn):
            n1 = nib.load(n1_fn) # preprocessed session 1
            n2 = nib.load(n2_fn) # preprocessed session 2
            
            bold1 = n1.get_data()
            bold2 = n2.get_data()
            bold = np.concatenate([bold1, bold2],axis=-1)
            
            nii_out = nib.Nifti1Image(bold, affine=n1.affine, header=n1.header) # pass affine and header from last MNI image
            nii_out.set_data_dtype(np.float32)
            nib.save(nii_out, mri_out)
        
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
        else:
            print('{} {} missing task-{}!'.format(self.subject, self.session, task))
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
            print('{} {} missing task-{}!'.format(self.subject, self.session, task))
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
            print('{} {} missing task-{}!'.format(self.subject, self.session, task))
           
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
    
    # def rsa_run_mask(self, task='rsa'):
    #     """Equalize number of runs per session per participant (in case of missing runs).
    #
    #     Args:
    #         task (str): which task. Default 'rsa'
    #
    #     Returns:
    #         rsa_runs (list): list of runs to use.
    #     """
    #     # check how many runs are in each session.
    #     ses1_mask = [
    #         os.path.exists(os.path.join(self.fmriprep_dir, '{}'.format(self.subject), '{}'.format('ses-mri01'), 'func', '{}_{}_task-cmrr2isomb4TR1500RSA_dir-AP_{}_space-T1w_desc-preproc_bold.nii.gz'.format(self.subject, 'ses-mri01', 'run-1'), )), # preprocessed run 1
    #
    #
    #         os.path.exists(os.path.join(self.fmriprep_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-01', task, 'run-02'), 'filtered_func_data.nii.gz')), # preprocessed run 2
    #         os.path.exists(os.path.join(self.fmriprep_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-01', task, 'run-03'), 'filtered_func_data.nii.gz')), # preprocessed run 3
    #         os.path.exists(os.path.join(self.fmriprep_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-01', task, 'run-04'), 'filtered_func_data.nii.gz')), # preprocessed run 4
    #     ]
    #
    #     ses2_mask = [
    #         os.path.exists(os.path.join(self.fmriprep_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-02', task, 'run-01'), 'filtered_func_data.nii.gz')), # preprocessed run 1
    #         os.path.exists(os.path.join(self.fmriprep_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-02', task, 'run-02'), 'filtered_func_data.nii.gz')), # preprocessed run 2
    #         os.path.exists(os.path.join(self.fmriprep_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-02', task, 'run-03'), 'filtered_func_data.nii.gz')), # preprocessed run 3
    #         os.path.exists(os.path.join(self.fmriprep_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-02', task, 'run-04'), 'filtered_func_data.nii.gz')), # preprocessed run 4
    #     ]
    #
    #     # delete after debugging!
    #     ses2_mask = [1, 1, 1, 1]
    #
    #     # use the same runs for each session for this participant.
    #     rsa_mask = np.array(ses1_mask)*np.array(ses2_mask)
    #     rsa_mask = np.array(rsa_mask, dtype=bool) # bools crucial for masking!
    #
    #     rsa_runs = np.array(['run-01', 'run-02', 'run-03', 'run-04'])
    #     rsa_runs = rsa_runs[rsa_mask]
    #     print('{} rsa runs = {}'.format(self.subject, rsa_runs))
    #     return rsa_runs
        
    
    def rsa_combine_epi(self, task='rsa'):
        """Concatenate the 4 runs per session of EPI data to perform a single GLM (RSA task).
        
        Args:
            task (str): which task. Default 'rsa'
        
        Notes:
            The output is the concantenated bold of all runs (input to first level).
            TO DO: Equalize number of runs per session per participant (in case of missing runs).
        """
                
        # rsa_runs = self.rsa_run_mask() # check which runs to use
        preprocessed_tag = 'space-T1w_desc-preproc_bold'
        
        for session in ['ses-mri01','ses-mri02']:
            
            df_trs = pd.DataFrame() # save number of TRs for events file concatenation
            
            # output is the concantenated bold of all runs per session (input to first level)
            out_file = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_run-concat_{}.nii.gz'.format(self.subject, session, task, preprocessed_tag))
            
            bold = []
            for r, this_run in enumerate(['run-1', 'run-2', 'run-3', 'run-4']): 
                for this_file in os.listdir(os.path.join(self.fmriprep_dir, self.subject, session, 'func')):
                    if (this_run in this_file) and ('space-T1w_desc-preproc_bold.nii.gz' in this_file):
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
        
        # rsa_runs = self.rsa_run_mask() # check which runs to use
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


    def rsa_nuisance_regressors(self, task='rsa'):
        """Concatenate the nifti files 4 runs of motion parameters, RETROICOR regressors, and 4th ventricle, from preprocessing.
        
        Args:
            task (str): which task. Default 'rsa'.
        
        Notes:
            Also makes the run regressors as single column EV text files.
            Equalize number of runs per session per participant (in case of missing runs).
        """

        rsa_runs = self.rsa_run_mask() # check which runs to use
        
        for self.session in ['ses-01','ses-02']:
            
            ####### MOTION PARAMETERS #######
            for mcf in np.arange(6): # 6 motion regressors
                # output is the concantenated niftis of all runs per session (nuisance input to first level)
                out_file = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_mcf_par{}.nii.gz'.format(self.subject, self.session, task, mcf + 1))
                
                mcf_nii = [] # all runs, current motion parameter
                for this_run in rsa_runs:
                    
                    nii = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'mc', 'mcf_par{}.nii.gz'.format(mcf + 1)))
                    mcf_nii.append(nii.get_data())

                mcf_nii = np.concatenate(mcf_nii, axis=-1)
            
                n1 = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'mc', 'mcf_par{}.nii.gz'.format(1)))
                out_data = nib.Nifti1Image(mcf_nii, affine=n1.affine, header=n1.header) # pass affine and header from last MNI image
                out_data.set_data_dtype(np.float32)
                nib.save(out_data, out_file)
                print(mcf_nii.shape)
            
            ####### RETROICOR #######
            for phys in np.arange(20): # 6 motion regressors
                # output is the concantenated niftis of all runs per session (nuisance input to first level)
                out_file = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_pnm_ev{}.nii.gz'.format(self.subject, self.session, task, phys + 1))
                
                phys_nii = [] # all runs, current motion parameter
                for this_run in rsa_runs:
                    nii = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'pnm_ev{:03}.nii.gz'.format(phys + 1)))
                    phys_nii.append(nii.get_data())

                phys_nii = np.concatenate(phys_nii, axis=-1)
            
                n1 = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'pnm_ev{:03}.nii.gz'.format(1)))
                out_data = nib.Nifti1Image(phys_nii, affine=n1.affine, header=n1.header) # pass affine and header from last MNI image
                out_data.set_data_dtype(np.float32)
                nib.save(out_data, out_file)
                print(phys_nii.shape)
            
            ####### 4th ventricle #######
            # output is the concantenated niftis of all runs per session (nuisance input to first level)
            out_file = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_pnm_ev4th_ventricle.nii.gz'.format(self.subject, self.session, task))
            
            phys_nii = [] # all runs, current motion parameter
            for this_run in rsa_runs:
                nii = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'pnm_ev4th_ventricle.nii.gz'))
                phys_nii.append(nii.get_data())

            phys_nii = np.concatenate(phys_nii, axis=-1)
        
            n1 = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'pnm_ev4th_ventricle.nii.gz'))
            out_data = nib.Nifti1Image(phys_nii, affine=n1.affine, header=n1.header) # pass affine and header from last MNI image
            out_data.set_data_dtype(np.float32)
            nib.save(out_data, out_file)
            print(phys_nii.shape)
            
            ####### Run regressors #######
            
            run_reg = []
            for this_run in rsa_runs:
                # just use 1st motion parameter to get number of volumes per run
                nii = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'mc', 'mcf_par{}.nii.gz'.format(mcf + 1)))
                run_reg.append(np.repeat(this_run, nii.shape[-1]))
            run_reg = np.concatenate(run_reg, axis=-1)
            df = pd.DataFrame(run_reg, columns=['run_regressors']) # make data frame to get_dummies
            df_dummies = df['run_regressors'].str.get_dummies() # splits single column into dummy variable columns for each individual run.
            
            for reg in df_dummies.columns.values:
                # text_file = open(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_confound_{}.txt'.format(self.subject, self.session, task, reg)), 'w')
                # df_dummies[reg].to_csv(text_file, header=None, index=None)
                # text_file.close()
                
                # save as a nifti file for regression
                reg_data = np.reshape(np.array(df_dummies[reg]), (1,1,1,run_reg.shape[-1]))
                nii_file = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_confound-{}.nii.gz'.format(self.subject, self.session, task, reg))
                reg_data = nib.Nifti1Image(reg_data, nii.affine, nii.header)
                reg_data.set_data_dtype(np.float32)
                nib.save(reg_data, nii_file)
                print(reg_data.shape)
                
            shell()

        print('success: rsa_nuisance_regressors {}'.format(self.subject))
        
        
    def nuisance_regressor_list(self, task):
        """Create a list of all confound file nuisance regressors for 1st level analysis.
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            RSA runs:
                Physiological noise regressors: 20
                4th ventricle: 1
                Motion regressors: 6
                Run regressors: 3-4 depending on missing data
        
            Localizers:
                Motion regressors: 6
        """
        
        if task in ['letters', 'colors']:
            
            feat_dir = os.path.join(self.first_level_dir, 'task-{}'.format(task))            
            text_file = open(os.path.join(feat_dir, '{}_task-{}_confound_evs_list.txt'.format(self.subject, task)), 'w')
            
            # grab motion regressors:
            regressors = [reg for reg in np.sort(glob.glob(os.path.join(feat_dir, '{}_task-{}_mcf_par*.nii.gz'.format(self.subject, task))))]
            for reg in regressors:
                text_file.write('{}\n'.format(reg))
            text_file.close()
            
        else: # rsa task
            
            for self.session in ['ses-01', 'ses-02']:
                
                feat_dir = os.path.join(self.first_level_dir, 'task-rsa')        
                text_file = open(os.path.join(feat_dir, '{}_{}_task-rsa_confound_evs_list.txt'.format(self.subject, self.session)), 'w')
        
                # grab RETROICOR regressors and 4th ventricle:
                regressors = [reg for reg in glob.glob(os.path.join(feat_dir, '{}_{}_task-rsa_pnm_ev*.nii.gz'.format(self.subject, self.session)))]
                for reg in regressors:
                    text_file.write('{}\n'.format(reg))
        
                # grab motion regressors:
                regressors = [reg for reg in np.sort(glob.glob(os.path.join(feat_dir, '{}_{}_task-rsa_mcf_par*.nii.gz'.format(self.subject, self.session))))]
                for reg in regressors:
                    text_file.write('{}\n'.format(reg))
                    
                # grab motion regressors:
                regressors = [reg for reg in glob.glob(os.path.join(feat_dir, '{}_{}_task-rsa_confound-run*.nii.gz'.format(self.subject, self.session)))]
                for reg in regressors:
                    text_file.write('{}\n'.format(reg))
                    
                text_file.close()        
        
        print('success: nuisance_regressor_list')


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

        for self.session in ['ses-01','ses-02']: # load session data with onsets
            events = pd.read_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_events.tsv'.format(self.subject, self.session, task)), sep='\t') # save concantenated events file
            
            #######################
            ### odd ball trials ###
            out_file = os.path.join(self.timing_files_dir, 'task-rsa', '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'oddballs'))
            # main regressors
            first = np.array(events[(events['oddball']==1)]['onset']) # onset in s
            second = events[(events['oddball']==1)]['RT'] # duration in s
            # replace any missed trials with stimulus duration
            idx = second.isnull() # doesn't recognize np.nan, only works on pd.series
            second[idx] = np.unique(events['duration'])[0] # stimulus duration seconds
            third = np.array(np.repeat(1, len(first)), dtype=int) # amplitude
            output = np.array(np.vstack((first, np.array(second), third)).T) # 1 x 3
            np.savetxt(out_file, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
            print(out_file)
            
            #######################
            #### drop oddballs ####
            events.drop(events[events['oddball'] == 1].index, inplace=True)

            # generate 3 column files for each of the 2x2 conditions
            for l,lcond in enumerate(np.unique(events['letter'])): # letters
                for c,ccond in enumerate(color_name): # 13 trained colors
                    # this letter in this color?
                    this_ev = events[(events['letter']==lcond) & (events['color_name']==ccond)]
                    # check if letter-color condition exists to prevent writing out empty EV timing files
                    if not this_ev.empty: 
                        # all files have to have the same names for subjects
                        out_file = os.path.join(self.timing_files_dir, 'task-rsa', '{}_{}_task-{}_{}_{}.txt'.format(self.subject, self.session, task, lcond, np.array(this_ev['trial_type_color'])[0]))
                        # main regressors
                        first = np.array(this_ev['onset']) # onset in s
                        second = np.array(this_ev['duration']) # onset in s # duration in s
                        third = np.array(np.repeat(1, len(first)),dtype=int) # amplitude
                        output = np.array(np.vstack((first, second, third)).T) # 1 x 3
                        np.savetxt(out_file, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
                        print(out_file)
        print('success: rsa_timing_files_letters')    
        
        
    def rsa_letters_fsf(self, task='rsa'):
        """Create the FSF files for each subject's first level analysis - RSA design each letter in each color.
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            Event-related design.
            Oddball trials are nuisance EV, physiological and motion regressors are listed in the confound text file.
            Run the actual FSF as batch script or from the command line: feat task-rsa_sub-01_ses-01.fsf
        """

        template_filename = os.path.join(self.analysis_dir, 'templates', 'task-{}_letters_first_level_template.fsf'.format(task))

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

        for self.session in ['ses-01','ses-02']:
            fsf_filename = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_letters.fsf'.format(self.subject, self.session, task)) # save fsf
            output_path = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_letters'.format(self.subject, self.session, task))

            bold = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_filtered_func_data.nii.gz'.format(self.subject, self.session, task))
            # calculate size of input data
            nii = nib.load(bold).get_data() # only do once
            nr_trs = str(nii.shape[-1])
            nr_voxels = str(nii.size)

            # motion parameters, retroicor, 4th ventricle, and columns for each run
            nuisance_regressors = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_confound_evs_list.txt'.format(self.subject, self.session, task))

            # timing files for each EV
            # black
            EV1_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'a_black'))
            EV2_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'b_black'))
            EV3_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'c_black'))
            EV4_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'d_black'))
            EV5_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'e_black'))
            EV6_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'f_black'))
            EV7_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'g_black'))
            EV8_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'h_black'))
            EV9_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'i_black'))
            EV10_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'j_black'))
            EV11_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'k_black'))
            EV12_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'l_black'))
            EV13_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'m_black'))
            EV14_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'n_black'))
            EV15_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'o_black'))
            EV16_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'p_black'))
            EV17_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'q_black'))
            EV18_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'r_black'))
            EV19_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'s_black'))
            EV20_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'t_black'))
            EV21_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'u_black'))
            EV22_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'v_black'))
            EV23_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'w_black'))
            EV24_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'x_black'))
            EV25_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'y_black'))
            EV26_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'z_black'))
            # color
            EV27_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'a_color'))
            EV28_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'b_color'))
            EV29_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'c_color'))
            EV30_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'d_color'))
            EV31_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'e_color'))
            EV32_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'f_color'))
            EV33_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'g_color'))
            EV34_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'h_color'))
            EV35_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'i_color'))
            EV36_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'j_color'))
            EV37_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'k_color'))
            EV38_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'l_color'))
            EV39_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'m_color'))
            EV40_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'n_color'))
            EV41_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'o_color'))
            EV42_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'p_color'))
            EV43_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'q_color'))
            EV44_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'r_color'))
            EV45_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'s_color'))
            EV46_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'t_color'))
            EV47_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'u_color'))
            EV48_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'v_color'))
            EV49_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'w_color'))
            EV50_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'x_color'))
            EV51_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'y_color'))
            EV52_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'z_color'))
            # oddballs last
            EV53_path = os.path.join(self.timing_files_dir,'task-{}'.format(task),'{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task,'oddballs'))

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
    
            # open preprocessing job and write command as new line
            cmd = 'feat {}'.format(fsf_filename)
            self.first_level_job = open(self.first_level_job_path, "a") # append is important, not write
            self.first_level_job.write(cmd)   # feat command
            self.first_level_job.write("\n\n")  # new line
            self.first_level_job.close()
        print('success: rsa_letters_fsf {}'.format(fsf_filename))


    def rsa_timing_files_2x2(self, task='rsa'):
        """Create GLM timing files for RSA - simple 2x2 comparison: trained/untrained vs. color/black
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            Event-related design.
            Outputs timing file for the oddball trials as well (nuisance EV).
        """
        
        for self.session in ['ses-01','ses-02']:
            events = pd.read_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task),'{}_{}_task-{}_events.tsv'.format(self.subject, self.session, task)),sep='\t') # save concantenated events file
            
            #######################
            ### odd ball trials ###
            out_file = os.path.join(self.timing_files_dir, 'task-rsa', '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'oddballs'))
            if os.path.exists(out_file):
                # concantenated file
                # main regressors
                first = np.array(events[(events['oddball']==1)]['onset']) # onset in s
                second = events[(events['oddball']==1)]['RT'] # duration in s
                # replace any missed trials with stimulus duration
                idx = second.isnull() # doesn't recognize np.nan, only works on pd.series
                second[idx] = np.unique(events['duration'])[0] # stimulus duration seconds
                third = np.array(np.repeat(1, len(first)),dtype=int) # amplitude
                output = np.array(np.vstack((first, np.array(second), third)).T) # 1 x 3
                np.savetxt(out_file, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
                print(out_file)
            
            # generate 3 column files for each of the 2x2 conditions
            for l,lcond in enumerate(['trained','untrained']): # letter condition
                for c,ccond in enumerate(['black','color']): # color condition
        
                    out_file = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_{}.txt'.format(self.subject, self.session, task, lcond, ccond))
                    # main regressors
                    first = np.array(events[(events['trial_type_letter']==lcond) & (events['trial_type_color']==ccond)]['onset']) # onset in s
                    second = np.array(events[(events['trial_type_letter']==lcond) & (events['trial_type_color']==ccond)]['duration']) # duration in s
                    third = np.array(np.repeat(1, len(first)),dtype=int) # amplitude
                    output = np.array(np.vstack((first, second, third)).T) # 1 x 3
                    np.savetxt(out_file, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
                    print(out_file)
        print('success: rsa_timing_files_2x2')    
        
        
    def rsa_2x2_fsf(self, task='rsa'):
        """Createf the FSF files for each subject's first level analysis - RSA 2x2 trained/untrained vs. black/color
                
        Args:
            task (str): which task. 
        
        Notes:
        ------
            Event-related design.
            Oddball trials are nuisance EV, physiological and motion regressors are listed in the confound text file.
            Run the actual FSF as batch script or from the command line: feat task-rsa_sub-01_ses-01.fsf
        """
        
        template_filename = os.path.join(self.analysis_dir, 'templates', 'task-rsa_2x2_first_level_template.fsf')
    
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
        
        for self.session in ['ses-01','ses-02']:
            fsf_filename = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_2x2.fsf'.format(self.subject, self.session, task)) # save fsf
            output_path = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_2x2'.format(self.subject, self.session, task)) 
        
            bold = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_filtered_func_data.nii.gz'.format(self.subject, self.session, task)) 
            # calculate size of input data
            nii = nib.load(bold).get_data() # only do once 
            nr_trs = str(nii.shape[-1])
            nr_voxels = str(nii.size)
        
            # motion parameters and columns for each session's mean
            nuisance_regressors = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_confound_evs_list.txt'.format(self.subject, self.session, task))
        
            # timing files for each EV
            EV1_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'untrained_black'))
            EV2_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'untrained_color'))
            EV3_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'trained_black'))
            EV4_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'trained_color'))
            EV5_path = os.path.join(self.timing_files_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}.txt'.format(self.subject, self.session, task, 'oddballs'))
                
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
                filedata = filedata.replace(this_string,replacements[st])

            # write output file
            f = open(fsf_filename,'w')
            f.write(filedata)
            f.close()
    
            # open preprocessing job and write command as new line
            cmd = 'feat {}'.format(fsf_filename)
            self.first_level_job = open(self.preprocessing_job_path, "a") # append is important, not write
            self.first_level_job.write(cmd)   # feat command
            self.first_level_job.write("\n\n")  # new line
            self.first_level_job.close()
        print('success: rsa_2x2_fsf {}'.format(fsf_filename))


