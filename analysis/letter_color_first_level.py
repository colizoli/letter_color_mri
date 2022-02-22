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
from IPython import embed as shell # for Oly's debugging only

class preprocess_class(object):
    def __init__(self, subject, session, analysis_dir, deriv_dir, mask_dir, template_dir, timing_files_dir):        
        self.subject        = str(subject)
        self.session        = session
        self.analysis_dir   = str(analysis_dir)
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.timing_files_dir = str(timing_files_dir)
        
        if not os.path.isdir(self.source_dir):
            os.mkdir(self.source_dir)
        if not os.path.isdir(self.deriv_dir):
            os.mkdir(self.deriv_dir)
        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
        if not os.path.isdir(self.timing_files_dir):
            os.mkdir(self.timing_files_dir)   
    

    def loc_timing_files(task, subject, session, runs):
        # GLM timing files for localizers
        # ABAB blocked designs
    
        stim_dur = 0.75 # stimulus duration seconds
    
        for run in runs:
            onsets = pd.read_csv(os.path.join(deriv_dir,subject,session,'events','{}_task-{}_{}_run-{}.csv'.format(subject,task, session,run)))
            onsets.drop(['Unnamed: 0'], axis=1, inplace=True) 
    
            for c,cond in enumerate(np.unique(onsets['condition'])):
                outFile = os.path.join(deriv_dir,'timing_files','task-{}_{}_{}_run-{}_{}.txt'.format(task,subject,session,run,cond))

                # main regressors
                first = np.array(onsets[onsets['condition']==cond]['onset_time'])
                second = np.repeat(stim_dur, len(first))
                third = np.array(np.repeat(1, len(first)),dtype=int)
                output = np.array(np.vstack((first, second, third)).T ,dtype=float) # 1 x 3
                f = open(outFile,'w')
                np.savetxt(f, output, fmt=['%.2f', '%.2f', '%i']) # doesn't work python 3
            
                f.close()
        
        print('success: loc_timing_files {} {} {} run-{}'.format(task,subject,session,run))     

    def rsa_timing_files(task, subject, session, runs):
        # GLM timing files for localizers
        # ABAB blocked designs
    
        stim_dur = 1.5 # stimulus duration seconds
    
        for run in runs:
            onsets = pd.read_csv(os.path.join(deriv_dir,subject,session,'events','{}_task-{}_{}_run-{}.csv'.format(subject,task, session,run)))
            onsets.drop(['Unnamed: 0'], axis=1, inplace=True) 
    
            for c,cond in enumerate(np.unique(onsets['color_condition'])):
                outFile = os.path.join(deriv_dir,'timing_files','task-{}_{}_{}_run-{}_{}.txt'.format(task,subject,session,run,cond))

                # main regressors
                first = np.array(onsets[onsets['color_condition']==cond]['onset_time'])
                second = np.repeat(stim_dur, len(first))
                third = np.repeat(1, len(first))
                output = np.vstack((first, second, third)).T # 1 x 3
                f = open(outFile,'w')
                np.savetxt(f, output, fmt=['%.2f', '%.2f', '%i'])
                f.close()
        
        print('success: rsa_timing_files {} {} {} run-{}'.format(task,subject,session,run))    

    def loc_fsf(task, subject, session, runs):
        # Creates the FSF files for each subject's first level analysis
        # localizers

        rm_TRs = str(5) # TRs to trim in beginning
        fwhm = str(3) # smoothing kernel mm
    
        # template_filename = os.path.join(deriv_dir,'first_level','templates','task-{}_template.fsf'.format(task))
        template_filename = os.path.join(deriv_dir,'first_level','templates','task-{}_template_B0_False.fsf'.format(task))
    
        markers = [
            '[$OUTPUT_PATH]', 
            '[$TR]',
            '[$NR_TRS]', 
            '[$REMOVE_TRS]',    # number of TRs to remove
            '[$DWELL_TIME]',    # EPI DWELL TIME, EFFECTIVE ECHO SPACING
            '[$ECHO_TIME]',     # EPI echo time
            '[$FWHM]',          # smoothing kernel
            '[$INPUT_FILENAME]', 
            '[$FMAP]',
            '[$FMAP_MAG_BRAIN]',
            '[$T1]',
            '[$EV1_FILENAME]',
            '[$EV2_FILENAME]', 
            '[$NR_VOXELS]'
        ]
    
        for r,run in enumerate(runs):        
            BOLD = os.path.join(deriv_dir,subject,session,'func','{}_task-{}_run-{}_bold.nii.gz'.format(subject,task,run))
            # calculate size of input data
            nii = nib.load(BOLD).get_data() # only do once 
            nr_trs = str(nii.shape[-1])
            nr_voxels = str(nii.size)
                
            FSF_filename = os.path.join(deriv_dir,'first_level','task-{}'.format(task),'task-{}_{}_{}_run-{}.fsf'.format(task,subject,session,run) ) # save fsf
            # replacements
            output_path = os.path.join(deriv_dir,'first_level','task-{}'.format(task),'task-{}_{}_{}_run-{}'.format(task,subject,session,run) ) 
            EV1_path = os.path.join(deriv_dir,'timing_files','task-{}_{}_{}_run-{}_{}.txt'.format(task,subject,session,run,'black'))
            EV2_path = os.path.join(deriv_dir,'timing_files','task-{}_{}_{}_run-{}_{}.txt'.format(task,subject,session,run,'color'))
    
            T1 = os.path.join(deriv_dir,subject,session,'anat','{}_T1w_brain.nii.gz'.format(subject))
            FMAP = os.path.join(deriv_dir,subject,session,'fmap','{}_acq-fmap.nii.gz'.format(subject))
            FMAP_MAG_BRAIN = os.path.join(deriv_dir,subject,session,'fmap','{}_fmap_run-01_fmap_brain.nii.gz'.format(subject))
            dwell_time = str(effective_EPI_echo_spacing[r])
            echo_time = str(echo_times[r])
            TR = str(TRs[r])
        
            replacements = [ # needs to match order of 'markers'
                output_path,
                TR,
                nr_trs,
                rm_TRs,
                dwell_time,
                echo_time,
                fwhm,
                BOLD,
                FMAP,
                FMAP_MAG_BRAIN,
                T1,
                EV1_path,
                EV2_path,
                nr_voxels
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
            print('success: {}_fsf {} {} run-{}'.format(task,subject,session,run))  
    
# ----------------------
# RUN FUNCTIONS
# ----------------------
subjects = ['sub-503']
sessions = ['sess-01']
color_runs = [
    ['01','02'] # multiband6, multiband4
]
rsa_runs = [
    ['1','2']
]
letter_runs = [
    ['1','2'] # multiband6, multiband4
]
vwfa_runs = [
    ['1','2'] # multiband6, multiband4
]

# TRs = [1.0,1.5]
TRs = [1.5,1.5]


#### B0 unwarping ####
# current issue: FEAT asks for echo time to run unwarping, but what is the echo time for the combined images?
# however, when unwarping run from FUGUE, it only asks for dwell time and magnitude image.
# what is FEAT doing differently from FUGUE? (combining with motion correct and BBR, but why is echo time needed?)
# 'dwell_time' refers to effective echo spacing of EPI data
echo_tdiff = 0.00246 # seconds, difference between echo times of field map magnitude images [0.0047,0.00716]
echo_times = [34.00, 39.60]
effective_EPI_echo_spacing = [0.000619987,0.000580009]*1000 # seconds, 'dwell time', param in json, corrected for phase encoding acceleration (not multiband acceleration)
unwarp_dir = 'y-' # A >> P, 'y-' as input
signal_loss_thresh = 10 # % default value


for s,subject in enumerate(subjects):
    for session in sessions:
        loc_timing_files('letters', subject, session, letter_runs[s]) # color localizer
        loc_timing_files('vwfa', subject, session, vwfa_runs[s]) # color localizer
        rsa_timing_files('rsa', subject, session, rsa_runs[s]) # rsa, run as localizer color vs. black for pilot
        loc_fsf('letters', subject, session, letter_runs[s])
        loc_fsf('vwfa', subject, session, vwfa_runs[s])
        loc_fsf('rsa', subject, session, rsa_runs[s])
        






