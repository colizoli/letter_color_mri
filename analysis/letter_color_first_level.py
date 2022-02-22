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
    def __init__(self, subject, session, analysis_dir, deriv_dir, mask_dir, template_dir, timing_files_dir, TR, FWHM, EPI_TE, EPI_EECHO):        
        self.subject        = 'sub-'+str(subject)
        self.session        = str(session)
        self.analysis_dir   = str(analysis_dir)
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.timing_files_dir = str(timing_files_dir)
        self.TR             = str(TR)
        self.FWHM           = str(FWHM)
        self.EPI_TE         = str(EPI_TE)
        self.EPI_EECHO      = str(EPI_EECHO)
        
        self.first_level_dir = os.path.join(self.deriv_dir,'first_level')
        
        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
        if not os.path.isdir(self.timing_files_dir):
            os.mkdir(self.timing_files_dir)   
            
        if not os.path.isdir(self.first_level_dir):
            os.mkdir(self.first_level)
            
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-colors')):
            os.mkdir(os.path.join(self.first_level_dir,'task-colors'))
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-letters')):
            os.mkdir(os.path.join(self.first_level_dir,'task-letters'))
        if not os.path.isdir(os.path.join(self.first_level_dir,'task-rsa')):
            os.mkdir(os.path.join(self.first_level_dir,'task-rsa'))
            
    def loc_timing_files(self, task):
        # GLM timing files for localizers
        # ABAB blocked designs
    
        stim_dur = 0.75 # stimulus duration seconds
    
        events = pd.read_csv(os.path.join(self.deriv_dir,self.subject,self.session,'func','{}_{}_task-{}loc_events.tsv'.format(self.subject,self.session,task)),sep='\t')
        
        for c,cond in enumerate(np.unique(events['trial_type'])):
            outFile = os.path.join(self.deriv_dir,'timing_files','task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,cond))

            # main regressors
            first = np.array(events[events['trial_type']==cond]['onset'])
            second = np.repeat(stim_dur, len(first))
            third = np.array(np.repeat(1, len(first)),dtype=int)
            output = np.array(np.vstack((first, second, third)).T) # 1 x 3
            np.savetxt(outFile, output, delimiter='/t', fmt='%.2f %.2f %i') #3rd column has to be an integer for FSL!
            print(outFile)

        print('success: loc_timing_files {} {} {}'.format(task,self.subject,self.session))     

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

    def loc_fsf(self,task):
        # Creates the FSF files for each subject's first level analysis - localizers
            
        template_filename = os.path.join(self.analysis_dir,'templates','task-{}loc_first_level_template.fsf'.format(task))
    
        markers = [
            '[$OUTPUT_PATH]', 
            '[$TR]',
            '[$NR_TRS]',        # number of volumes
            '[$NR_VOXELS]',     # total number of voxels
            '[$EPI_EECHO]',     # EPI DWELL TIME, EFFECTIVE ECHO SPACING
            '[$EPI_TE]',        # EPI echo time
            '[$FWHM]',          # smoothing kernel
            '[$INPUT_FILENAME]', 
            '[$FMAP]',
            '[$FMAP_MAG_BRAIN]',
            '[$T1_BRAIN]',
            '[$EV1_FILENAME]',
            '[$EV2_FILENAME]', 
        ]
    
        BOLD = os.path.join(self.deriv_dir,self.subject,self.session,'func','{}_{}_task-{}_bold.nii.gz'.format(self.subject,self.session,task))
        # calculate size of input data
        nii = nib.load(BOLD).get_data() # only do once 
        nr_trs = str(nii.shape[-1])
        nr_voxels = str(nii.size)
            
        FSF_filename = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_{}.fsf'.format(task,self.subject,self.session) ) # save fsf
        # replacements
        output_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_{}_{}'.format(task,self.subject,self.session)) 
        EV1_path = os.path.join(self.deriv_dir,'timing_files','task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'Color'))
        EV2_path = os.path.join(self.deriv_dir,'timing_files','task-{}_{}_{}_{}.txt'.format(task,self.subject,self.session,'Black'))

        T1_BRAIN = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w_brain.nii.gz'.format(self.subject,self.session))
        FMAP = os.path.join(self.deriv_dir,self.subject,self.session,'fmap','{}_{}_acq-fmap.nii.gz'.format(self.subject,self.session))
        FMAP_MAG_BRAIN = os.path.join(self.deriv_dir,self.subject,self.session,'fmap','{}_{}_run-01_fmap_brain.nii.gz'.format(self.subject, self.session))
        
        replacements = [ # needs to match order of 'markers'
            output_path,
            self.TR,
            nr_trs,
            nr_voxels,
            self.EPI_EECHO, # dwell time is effective echo spacing (EPI data not field map!!)
            self.EPI_TE,
            self.FWHM,
            BOLD,
            FMAP,
            FMAP_MAG_BRAIN,
            T1_BRAIN,
            EV1_path,
            EV2_path,
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







