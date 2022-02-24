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

        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
        if not os.path.isdir(self.timing_files_dir):
            os.mkdir(self.timing_files_dir)   
        
        self.preprocess_dir = os.path.join(self.deriv_dir,'preprocessing')
        
        self.first_level_dir = os.path.join(self.deriv_dir,'first_level')
        if not os.path.isdir(self.first_level_dir):
            os.mkdir(self.first_level)
            
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-colors')):
            os.mkdir(os.path.join(self.first_level_dir,'task-colors'))
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-letters')):
            os.mkdir(os.path.join(self.first_level_dir,'task-letters'))
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-rsa')):
            os.mkdir(os.path.join(self.first_level_dir,'task-rsa'))
            
    def loc_combine_epi(self, task):
        # concatenate the 2 blocks of EPI data to perform a single GLM
        
        # output is the concantenated bold of both sessions (input to first level)
        outFile = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_bold_mni.nii.gz'.format(task,self.subject)) 
        N1 = nib.load(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,'ses-01'),'filtered_func_data_mni.nii.gz')) # preprocessed session 1
        N2 = nib.load(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,'ses-02'),'filtered_func_data_mni.nii.gz')) # preprocessed session 2
        BOLD1 = N1.get_data()
        BOLD2 = N2.get_data()
        
        BOLD = np.concatenate([BOLD1,BOLD2],axis=-1)
        outData = nib.Nifti1Image(BOLD, affine=N1.affine, header=N1.header) # pass affine and header from last MNI image
        outData.set_data_dtype(np.float32)
        nib.save(outData, outFile)
        
        print('success: concatenate_epi_data {}'.format(self.subject))
     
    def loc_combine_timing_files(self, task):
        # concatenate the timing files: 2nd session have to add time = #TRs * TR
        # GLM timing files for blocked design (localizers)
        # ABAB blocked designs
        stim_dur = 0.75 # stimulus duration seconds
        
        # take FIRST session's BOLD to count TRs to add to 2nd sessions' onsets
        BOLD1 = nib.load(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,'ses-01'),'filtered_func_data_mni.nii.gz'))
        ntrs = BOLD1.shape[-1]  # number of TRs
        time2add = ntrs*float(self.TR) # time to add in seconds to block 2's onsets
        print('Time to add to run 2: {}'.format(time2add))
        
        # open block 2's events
        events2 = pd.read_csv(os.path.join(self.deriv_dir,self.subject,'ses-02','func','{}_{}_task-{}_events.tsv'.format(self.subject,'ses-02',task)),sep='\t')
        events2['onset'] =  events2['onset'] + time2add
        # open block 1's events and concantenate
        events1 = pd.read_csv(os.path.join(self.deriv_dir,self.subject,'ses-01','func','{}_{}_task-{}_events.tsv'.format(self.subject,'ses-01',task)),sep='\t')
        events = pd.concat([events1,events2],axis=0)

        # generate 3 column files for each of the 2x2 conditions
        for c,cond in enumerate(np.unique(events['trial_type'])):
            outFile = os.path.join(self.deriv_dir,'timing_files','task-{}_{}_{}.txt'.format(task,self.subject,cond))

            # main regressors
            first = np.array(events[events['trial_type']==cond]['onset'])
            second = np.repeat(stim_dur, len(first))
            third = np.array(np.repeat(1, len(first)),dtype=int)
            output = np.array(np.vstack((first, second, third)).T) # 1 x 3
            np.savetxt(outFile, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
            print(outFile)
        print('success: combine_timing_files {}'.format(self.subject))
    
    def loc_nuisance_regressors(self, task):
        # concatenate the 2 sessions of motion parameters from preprocessing
        # these are found in derivatives/preprocessing/task/task_subject_session.feat/mc/prefiltered_func_data_mcf.par
        # Nrows = NTRs, Ncols = 6 (mc directions), note space separated
        # This function also outputs the columns of 1s and 0s for each blocks' mean 
        
        #### Motion parameters ####
        mc1 = pd.read_csv(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,'ses-01'),'mc','prefiltered_func_data_mcf.par'),header=None,sep='\s+',float_precision='round_trip')
        mc2 = pd.read_csv(os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,'ses-02'),'mc','prefiltered_func_data_mcf.par'),header=None,sep='\s+',float_precision='round_trip')
            
        for col in mc1.columns.values: # convert data to numeric
            mc1[col] = pd.to_numeric(mc1[col])
            mc2[col] = pd.to_numeric(mc2[col])

        mc = pd.concat([mc1,mc2],axis=0) # concantenate the motion regressors
        
        #### Session means - make columns of 1s and 0s for the length of each session ####
        b1 = np.concatenate((np.repeat(1,len(mc1)),np.repeat(0,len(mc2))),axis=0) # session 1: 1s then 0s
        b2 = np.concatenate((np.repeat(0,len(mc1)),np.repeat(1,len(mc2))),axis=0) # sessoin 2: 0s then 1s
        # add to motion dataframe
        mc['b1'] = b1
        mc['b2'] = b2
        
        # save without header or index! suppress scientific notation!
        mc.to_csv(os.path.join(self.timing_files_dir,'task-{}_{}_nuisance_regressors.txt'.format(task,self.subject)),header=None,index=False,sep=',',float_format='%.15f')        
        print('success: loc_nuisance_regressors {}'.format(self.subject))
        
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
        nuisance_regressors = os.path.join(self.timing_files_dir,'task-{}_{}_nuisance_regressors.txt'.format(task,self.subject))
        
        # timing files for each EV
        EV1_path = os.path.join(self.deriv_dir,'timing_files','task-{}_{}_{}.txt'.format(task,self.subject,'Color'))
        EV2_path = os.path.join(self.deriv_dir,'timing_files','task-{}_{}_{}.txt'.format(task,self.subject,'Black'))
        
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
        print('success: loc_fsf {}'.format(FSF_filename))
        
    
    
    
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









