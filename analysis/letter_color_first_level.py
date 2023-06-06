#!/usr/bin/env python
# encoding: utf-8

"""
first_level.py
Created by O.Colizoli, June-2019
Last update: 11-06-2019
Python version 2.7

The following packages need to be installed because they are called from the command line, but not imported:
fsl

"""

# TO DO:
# bids output files: sub, session, task, run, filetypebids coverted does latter...
# run-01 or run-1?? bids coverted does latter...

import os, subprocess, sys
import shutil as sh
import nibabel as nib
import pandas as pd
import numpy as np
from IPython import embed as shell # for Oly's debugging only

class first_level_class(object):
    def __init__(self, subject, analysis_dir, deriv_dir, mask_dir, template_dir, timing_files_dir, TR):        
        self.subject        = 'sub-'+str(subject)
        self.analysis_dir   = str(analysis_dir)
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.timing_files_dir = str(timing_files_dir)
        self.TR             = str(TR)

        if not os.path.isdir(self.timing_files_dir):
            os.mkdir(self.timing_files_dir)   
            os.mkdir(os.path.join(self.timing_files_dir,'task-colors'))
            os.mkdir(os.path.join(self.timing_files_dir,'task-letters'))
            os.mkdir(os.path.join(self.timing_files_dir,'task-rsa'))
            
        self.preprocess_dir = os.path.join(self.deriv_dir,'preprocessing')
        
        self.first_level_dir = os.path.join(self.deriv_dir,'first_level')
        if not os.path.isdir(self.first_level_dir):
            os.mkdir(self.first_level_dir)
            
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-colors')):
            os.mkdir(os.path.join(self.first_level_dir,'task-colors'))
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-letters')):
            os.mkdir(os.path.join(self.first_level_dir,'task-letters'))
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-rsa')):
            os.mkdir(os.path.join(self.first_level_dir,'task-rsa'))
        
        # write unix commands to job to run in parallel
        self.first_level_job_path = os.path.join(self.analysis_dir,'jobs','job_first_level_{}.txt'.format(self.subject))
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
    
    def rsa_run_mask(self, task='rsa'):
        """Equalize number of runs per session per participant (in case of missing runs).
        
        Args:
            task (str): which task. Default 'rsa'
        
        Returns:
            rsa_runs (list): list of runs to use.
        """
        # check how many runs are in each session.
        ses1_mask = [
            os.path.exists(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-01', task, 'run-01'), 'filtered_func_data.nii.gz')), # preprocessed run 1
            os.path.exists(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-01', task, 'run-02'), 'filtered_func_data.nii.gz')), # preprocessed run 2
            os.path.exists(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-01', task, 'run-03'), 'filtered_func_data.nii.gz')), # preprocessed run 3
            os.path.exists(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-01', task, 'run-04'), 'filtered_func_data.nii.gz')), # preprocessed run 4
        ]
        
        ses2_mask = [
            os.path.exists(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-02', task, 'run-01'), 'filtered_func_data.nii.gz')), # preprocessed run 1
            os.path.exists(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-02', task, 'run-02'), 'filtered_func_data.nii.gz')), # preprocessed run 2
            os.path.exists(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-02', task, 'run-03'), 'filtered_func_data.nii.gz')), # preprocessed run 3
            os.path.exists(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-02', task, 'run-04'), 'filtered_func_data.nii.gz')), # preprocessed run 4
        ]
        
        # delete after debugging!        
        ses2_mask = [1, 1, 1, 1] 
        
        # use the same runs for each session for this participant.
        rsa_mask = np.array(ses1_mask)*np.array(ses2_mask)
        rsa_mask = np.array(rsa_mask, dtype=bool) # bools crucial for masking!
        
        rsa_runs = np.array(['run-01', 'run-02', 'run-03', 'run-04'])
        rsa_runs = rsa_runs[rsa_mask]
        print('{} rsa runs = {}'.format(self.subject, rsa_runs))
        return rsa_runs
        
    
    def rsa_combine_epi(self, task='rsa'):
        """Concatenate the 4 runs per session of EPI data to perform a single GLM (RSA task).
        
        Args:
            task (str): which task. Default 'rsa'
        
        Notes:
            The output is the concantenated bold of all runs (input to first level).
            Also makes the run regressors as EV text files.
            Equalize number of runs per session per participant (in case of missing runs).
        """
        
        rsa_runs = self.rsa_run_mask() # check which runs to use
        
        for self.session in ['ses-01','ses-02']:
            # output is the concantenated bold of all runs per session (input to first level)
            out_file = os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_filtered_func_data.nii.gz'.format(self.subject, self.session, task))
            
            bold = []
            for this_run in rsa_runs:
                nii = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'filtered_func_data.nii.gz'))
                bold.append(nii.get_data())
                        
            bold = np.concatenate(bold, axis=-1)
            
            n1 = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, 'ses-01', task, 'run-01'), 'filtered_func_data.nii.gz'))
            out_data = nib.Nifti1Image(bold, affine=n1.affine, header=n1.header) # pass affine and header from last MNI image
            out_data.set_data_dtype(np.float32)
            nib.save(out_data, out_file)
            print(bold.shape)
            shell()
        print('success: rsa_combine_epi')
        
        
    def rsa_combine_events(self, task='rsa'):
        """For the RSA task, concantenate the events files of all runs and output in first_level directory.
        
        Args:
            task (str): which task. Default 'rsa'
        
        Notes:
            Equalize number of runs per session per participant (in case of missing runs).
        """
        
        rsa_runs = self.rsa_run_mask() # check which runs to use
        
        for self.session in ['ses-01','ses-02']:
            
            nruns = len(rsa_runs)
            
            time2add = []
            
            # check how long each run was and calculate seconds
            for r, this_run in enumerate(rsa_runs):
                
                bold = nib.load(os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, self.session, task, this_run), 'filtered_func_data.nii.gz'))
                ntrs = bold.shape[-1]  # number of TRs
                time2add.append(ntrs*float(self.TR)) # time to add in seconds to next run's onsets
                print('Time to add to run {}: {}'.format(r + 1, ntrs*float(self.TR)))
                
            time2add = np.cumsum(time2add) # cummulative sum per run
            
            # open events files and add times
            all_events = pd.DataFrame()
            for r, this_run in enumerate(rsa_runs):
                if r==0: # first run
                    # open file but don't add time
                    events = pd.read_csv(os.path.join(self.deriv_dir, self.subject, self.session, 'func', '{}_{}_task-{}_{}_events.tsv'.format(self.subject, self.session, task, this_run)), sep='\t')
                    all_events = pd.concat([all_events, events], axis=0)
                else:
                    events = pd.read_csv(os.path.join(self.deriv_dir, self.subject, self.session, 'func', '{}_{}_task-{}_{}_events.tsv'.format(self.subject, self.session, task, this_run)), sep='\t')
                    events['onset'] =  events['onset'] + time2add[r]
                    all_events = pd.concat([all_events, events], axis=0)
                print(all_events.shape)

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
            all_events.to_csv(os.path.join(self.first_level_dir, 'task-{}'.format(task), '{}_{}_task-{}_events.tsv'.format(self.subject, self.session, task)), sep='\t') 

        print('success: rsa_combine_events')    


    def rsa_nuisance_regressors(self,task='rsa'):
        # concatenate the 4 runs of motion parameters from preprocessing
        # these are found in derivatives/preprocessing/task/task_subject_session.feat/mc/prefiltered_func_data_mcf.par
        # Nrows = NTRs, Ncols = 6 (mc directions), note space separated
        # This function also outputs the columns of 1s and 0s for each blocks' mean 
        # Also outputs the odd ball button presses in a seperate 3 column format file (so FEAT will convolve with same HRF as used in main analysis)
        # RT on oddball trials set as duration (not amplitude), if no button press then set to stimulus duration
        
        for self.session in ['ses-01','ses-02']:
            
            phys_fn = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, self.session, task), 'pnm_ev{}.nii.gz'.format(phys + 1))
            mcf_fn = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_preprocessing.feat'.format(self.subject, self.session, task), 'mc', 'mcf_par{}.nii.gz'.format(mcf + 1))
            
            
            #### Motion parameters of each run's preprocessing ####
            mc1 = pd.read_csv(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}_{}.feat'.format(task,self.subject,self.session,'run-01'),'mc','prefiltered_func_data_mcf.par'),header=None,sep='\s+',float_precision='round_trip')
            mc2 = pd.read_csv(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}_{}.feat'.format(task,self.subject,self.session,'run-02'),'mc','prefiltered_func_data_mcf.par'),header=None,sep='\s+',float_precision='round_trip')
            mc3 = pd.read_csv(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}_{}.feat'.format(task,self.subject,self.session,'run-03'),'mc','prefiltered_func_data_mcf.par'),header=None,sep='\s+',float_precision='round_trip')
            mc4 = pd.read_csv(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}_{}.feat'.format(task,self.subject,self.session,'run-04'),'mc','prefiltered_func_data_mcf.par'),header=None,sep='\s+',float_precision='round_trip')
        
            for col in mc1.columns.values: # convert data to numeric
                mc1[col] = pd.to_numeric(mc1[col])
                mc2[col] = pd.to_numeric(mc2[col])
                mc3[col] = pd.to_numeric(mc3[col])
                mc4[col] = pd.to_numeric(mc4[col])
            
            mc = pd.concat([mc1,mc2,mc3,mc4],axis=0) # concantenate the motion regressors
        
            #### Run means - make columns of 1s and 0s for the length of each run ####
            b1 = np.concatenate( (np.repeat(1,len(mc1))   ,  np.repeat(0,len(mc2))   ,  np.repeat(0,len(mc3))  , np.repeat(0,len(mc4)) ),axis=0) # session 1: 1s run 1
            b2 = np.concatenate( (np.repeat(0,len(mc1))   ,  np.repeat(1,len(mc2))   ,  np.repeat(0,len(mc3))  , np.repeat(0,len(mc4)) ),axis=0) # session 2: 1s run 2
            b3 = np.concatenate( (np.repeat(0,len(mc1))   ,  np.repeat(0,len(mc2))   ,  np.repeat(1,len(mc3))  , np.repeat(0,len(mc4)) ),axis=0) # session 3: 1s run 3
            b4 = np.concatenate( (np.repeat(0,len(mc1))   ,  np.repeat(0,len(mc2))   ,  np.repeat(0,len(mc3))  , np.repeat(1,len(mc4)) ),axis=0) # session 4: 1s run 4
        
            # add to motion dataframe
            mc['b1'] = b1
            mc['b2'] = b2
            mc['b3'] = b3
            mc['b4'] = b4
        
            # save without header or index! suppress scientific notation!
            mc.to_csv(os.path.join(self.timing_files_dir,'task-{}'.format(task),'task-{}_{}_{}_nuisance_regressors.txt'.format(task,self.subject,self.session)),header=None,index=False,sep=',',float_format='%.15f')  
            #######################
            ### odd ball trials ###
            stim_dur = 1.5 # stimulus duration seconds
            # concantenated file
            events = pd.read_csv(os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_{}_events.tsv'.format(task,self.subject,self.session)),sep='\t') # save concantenated events file
            outFile = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'oddballs'))
            # main regressors
            first = np.array(events[(events['oddball']==1)]['onset']) # onset in s
            second = events[(events['oddball']==1)]['RT'] # duration in s
            # replace any missed trials with stimulus duration
            idx=second.isnull() # doesn't recognize np.nan, only works on pd.series
            second[idx] = stim_dur
            third = np.array(np.repeat(1, len(first)),dtype=int) # amplitude
            output = np.array(np.vstack((first, np.array(second), third)).T) # 1 x 3
            np.savetxt(outFile, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
            print(outFile)
        print('success: loc_nuisance_regressors {}'.format(self.subject))


    def rsa_timing_files_letters(self,task='rsa'):
        # GLM timing files for RSA - each letter in each color (2 trials per run)
        # event related design
        
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
            events = pd.read_csv(os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_{}_events.tsv'.format(task,self.subject,self.session)),sep='\t') 
            # drop oddballs        
            events.drop(events[events['oddball'] == 1].index, inplace=True)

            # generate 3 column files for each of the 2x2 conditions
            for l,lcond in enumerate(np.unique(events['letter'])): # letters
                for c,ccond in enumerate(color_name): # 13 trained colors
                    # this letter in this color?
                    this_ev = events[(events['letter']==lcond) & (events['color_name']==ccond)]
                    # check if letter-color condition exists to prevent writing out empty EV timing files
                    if not this_ev.empty: 
                        # all files have to have the same names for subjects
                        outFile = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}_{}.txt'.format(task,self.subject,self.session,lcond,np.array(this_ev['trial_type_color'])[0]))
                        # main regressors
                        first = np.array(this_ev['onset']) # onset in s
                        second = np.array(this_ev['duration']) # onset in s # duration in s
                        third = np.array(np.repeat(1, len(first)),dtype=int) # amplitude
                        output = np.array(np.vstack((first, second, third)).T) # 1 x 3
                        np.savetxt(outFile, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
                        print(outFile)
        print('success: rsa_timing_files_letters')    
        
    def rsa_letters_fsf(self,task='rsa'):
        # Creates the FSF files for each subject's first level analysis - RSA design each letter in each color
        # Run the actual FSF from the command line: feat task-rsa_sub-01_ses-01.fsf
            
        template_filename = os.path.join(self.analysis_dir,'templates','task-{}_letters_first_level_template.fsf'.format(task))

        markers = [
            '[$OUTPUT_PATH]',
            '[$NR_TRS]',
            '[$NR_VOXELS]',
            '[$INPUT_FILENAME]',
            '[$NUISANCE]',
            '[$MNI_BRAIN]',
            '[$EV1_FILENAME]','[$EV2_FILENAME]','[$EV3_FILENAME]','[$EV4_FILENAME]','[$EV5_FILENAME]','[$EV6_FILENAME]','[$EV7_FILENAME]','[$EV8_FILENAME]','[$EV9_FILENAME]', '[$EV10_FILENAME]',
            '[$EV11_FILENAME]','[$EV12_FILENAME]','[$EV13_FILENAME]','[$EV14_FILENAME]','[$EV15_FILENAME]','[$EV16_FILENAME]','[$EV17_FILENAME]','[$EV18_FILENAME]','[$EV19_FILENAME]', '[$EV20_FILENAME]',
            '[$EV21_FILENAME]','[$EV22_FILENAME]','[$EV23_FILENAME]','[$EV24_FILENAME]','[$EV25_FILENAME]','[$EV26_FILENAME]','[$EV27_FILENAME]','[$EV28_FILENAME]','[$EV29_FILENAME]', '[$EV30_FILENAME]',
            '[$EV31_FILENAME]','[$EV32_FILENAME]','[$EV33_FILENAME]','[$EV34_FILENAME]','[$EV35_FILENAME]','[$EV36_FILENAME]','[$EV37_FILENAME]','[$EV38_FILENAME]','[$EV39_FILENAME]', '[$EV40_FILENAME]',
            '[$EV41_FILENAME]','[$EV42_FILENAME]','[$EV43_FILENAME]','[$EV44_FILENAME]','[$EV45_FILENAME]','[$EV46_FILENAME]','[$EV47_FILENAME]','[$EV48_FILENAME]','[$EV49_FILENAME]', '[$EV50_FILENAME]',
            '[$EV51_FILENAME]','[$EV52_FILENAME]',
            '[$EV53_FILENAME]' # oddballs
        ]

        for self.session in ['ses-01','ses-02']:
            FSF_filename = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_letters_{}_{}.fsf'.format(task,self.subject,self.session)) # save fsf
            output_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_letters_{}_{}'.format(task,self.subject,self.session))

            BOLD = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_{}_bold_mni.nii.gz'.format(task,self.subject,self.session))
            # calculate size of input data
            nii = nib.load(BOLD).get_data() # only do once
            nr_trs = str(nii.shape[-1])
            nr_voxels = str(nii.size)

            # motion parameters and columns for each session's mean
            nuisance_regressors = os.path.join(self.timing_files_dir,'task-{}'.format(task),'task-{}_{}_{}_nuisance_regressors.txt'.format(task,self.subject,self.session))
            MNI_BRAIN  = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz')

            # timing files for each EV
            # black
            EV1_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'a_black'))
            EV2_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'b_black'))
            EV3_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'c_black'))
            EV4_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'d_black'))
            EV5_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'e_black'))
            EV6_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'f_black'))
            EV7_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'g_black'))
            EV8_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'h_black'))
            EV9_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'i_black'))
            EV10_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'j_black'))
            EV11_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'k_black'))
            EV12_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'l_black'))
            EV13_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'m_black'))
            EV14_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'n_black'))
            EV15_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'o_black'))
            EV16_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'p_black'))
            EV17_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'q_black'))
            EV18_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'r_black'))
            EV19_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'s_black'))
            EV20_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'t_black'))
            EV21_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'u_black'))
            EV22_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'v_black'))
            EV23_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'w_black'))
            EV24_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'x_black'))
            EV25_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'y_black'))
            EV26_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'z_black'))
            # color
            EV27_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'a_color'))
            EV28_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'b_color'))
            EV29_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'c_color'))
            EV30_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'d_color'))
            EV31_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'e_color'))
            EV32_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'f_color'))
            EV33_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'g_color'))
            EV34_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'h_color'))
            EV35_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'i_color'))
            EV36_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'j_color'))
            EV37_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'k_color'))
            EV38_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'l_color'))
            EV39_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'m_color'))
            EV40_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'n_color'))
            EV41_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'o_color'))
            EV42_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'p_color'))
            EV43_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'q_color'))
            EV44_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'r_color'))
            EV45_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'s_color'))
            EV46_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'t_color'))
            EV47_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'u_color'))
            EV48_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'v_color'))
            EV49_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'w_color'))
            EV50_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'x_color'))
            EV51_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'y_color'))
            EV52_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'z_color'))
            # oddballs last
            EV53_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'oddballs'))

            # replacements
            replacements = [ # needs to match order of 'markers'
                output_path,
                nr_trs,
                nr_voxels,
                BOLD,
                nuisance_regressors,
                MNI_BRAIN,
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
                filedata = filedata.replace(this_string,replacements[st])

            # write output file
            f = open(FSF_filename,'w')
            f.write(filedata)
            f.close()
    
            # open preprocessing job and write command as new line
            cmd = 'feat {}'.format(FSF_filename)
            self.first_level_job = open(self.first_level_job_path, "a") # append is important, not write
            self.first_level_job.write(cmd)   # feat command
            self.first_level_job.write("\n\n")  # new line
            self.first_level_job.close()
        print('success: rsa_letters_fsf {}'.format(FSF_filename))

    def rsa_timing_files_2x2(self,task='rsa'):
        # GLM timing files for RSA - simple 2x2 comparison: trained/untrained vs. color/black
        # event related design
        
        for self.session in ['ses-01','ses-02']:
            events = pd.read_csv(os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_{}_events.tsv'.format(task,self.subject,self.session)),sep='\t') # save concantenated events file
        
            # generate 3 column files for each of the 2x2 conditions
            for l,lcond in enumerate(['trained','untrained']): # letter condition
                for c,ccond in enumerate(['black','color']): # color condition
        
                    outFile = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}_{}.txt'.format(task,self.subject,self.session,lcond,ccond))
                    # main regressors
                    first = np.array(events[(events['trial_type_letter']==lcond) & (events['trial_type_color']==ccond)]['onset']) # onset in s
                    second = np.array(events[(events['trial_type_letter']==lcond) & (events['trial_type_color']==ccond)]['duration']) # duration in s
                    third = np.array(np.repeat(1, len(first)),dtype=int) # amplitude
                    output = np.array(np.vstack((first, second, third)).T) # 1 x 3
                    np.savetxt(outFile, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
                    print(outFile)
        print('success: rsa_timing_files_2x2')    
        
        
    def rsa_2x2_fsf(self,task='rsa'):
        # Creates the FSF files for each subject's first level analysis - RSA 2x2 trained/untrained vs. black/color
        # Run the actual FSF from the command line: feat task-rsa_sub-01_ses-01.fsf
            
        template_filename = os.path.join(self.analysis_dir,'templates','task-{}_2x2_first_level_template.fsf'.format(task))
    
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
            '[$MNI_BRAIN]'
        ]
        
        for self.session in ['ses-01','ses-02']:
            FSF_filename = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_2x2_{}_{}.fsf'.format(task,self.subject,self.session)) # save fsf
            output_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_2x2_{}_{}'.format(task,self.subject,self.session)) 
        
            BOLD = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_{}_bold_mni.nii.gz'.format(task,self.subject,self.session)) 
            # calculate size of input data
            nii = nib.load(BOLD).get_data() # only do once 
            nr_trs = str(nii.shape[-1])
            nr_voxels = str(nii.size)
        
            # motion parameters and columns for each session's mean
            nuisance_regressors = os.path.join(self.timing_files_dir,'task-{}'.format(task),'task-{}_{}_{}_nuisance_regressors.txt'.format(task,self.subject,self.session))
        
            # timing files for each EV
            EV1_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'untrained_black'))
            EV2_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'untrained_color'))
            EV3_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'trained_black'))
            EV4_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'trained_color'))
            EV5_path = os.path.join(self.deriv_dir,'timing_files','task-{}'.format(task),'task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'oddballs'))
        
            MNI_BRAIN  = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz')
        
            # replacements
            replacements = [ # needs to match order of 'markers'
                output_path,
                nr_trs,
                BOLD,
                nuisance_regressors,
                EV1_path,
                EV2_path,
                EV3_path,
                EV4_path,
                EV5_path,
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
    
            # open preprocessing job and write command as new line
            cmd = 'feat {}'.format(FSF_filename)
            self.first_level_job = open(self.preprocessing_job_path, "a") # append is important, not write
            self.first_level_job.write(cmd)   # feat command
            self.first_level_job.write("\n\n")  # new line
            self.first_level_job.close()
        print('success: rsa_2x2_fsf {}'.format(FSF_filename))


