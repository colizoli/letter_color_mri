#!/usr/bin/env python
# encoding: utf-8

"""
higher_level.py
Created by O.Colizoli, June-2019
Last update: 11-06-2019
Python version 2.7

The following packages need to be installed because they are called from the command line, but not imported:
fsl

"""

import os, subprocess, sys
import shutil as sh
import nibabel as nib
import pandas as pd
import numpy as np
from IPython import embed as shell # for Oly's debugging only

class higher_level_class(object):
    def __init__(self, subjects, sessions, analysis_dir, deriv_dir, mask_dir, template_dir, TR):        
        self.subjects       = subjects
        self.sessions       = sessions
        self.analysis_dir   = str(analysis_dir)
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.TR             = str(TR)

        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
            
        self.preprocess_dir = os.path.join(self.deriv_dir,'preprocessing')
        self.first_level_dir = os.path.join(self.deriv_dir,'first_level')
        
        self.higher_level_dir = os.path.join(self.deriv_dir,'higher_level')
        if not os.path.isdir(self.higher_level_dir):
            os.mkdir(self.higher_level_dir)
            
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-colors')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-colors'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-letters')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-letters'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-rsa')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-rsa'))
        
        # brain atlases 
        # note: 'maxprob-thr0-2mm' files have labels in single volume, while the 'prob-2mm' files contain probabilities for each label as time axis
        self.mni        = os.path.join(self.mask_dir,'MNI152_T1_2mm.nii.gz')
        self.mni_brain  = os.path.join(self.mask_dir,'MNI152_T1_2mm_brain.nii.gz')
        self.mni_labels = os.path.join(self.mask_dir,'atlases','MNI','MNI-maxprob-thr0-2mm.nii.gz')
        self.ho_cortical = os.path.join(self.mask_dir,'atlases','HarvardOxford','HarvardOxford-cort-maxprob-thr0-2mm.nii.gz')
        
        # write unix commands to job to run in parallel
        self.higher_level_job_path = os.path.join(self.analysis_dir,'jobs','job_higher_level.txt')
        if not os.path.exists(self.higher_level_job_path):
            self.higher_level_job_path = open(self.higher_level_job_path, "w")
            self.higher_level_job_path.write("#!/bin/bash\n")
            self.higher_level_job_path.close()
            
    def labels_harvard_oxford_whole_brain(self,):
        # combines the cortical and subcortical labels together in a single mask
        self.ho_cortical = os.path.join(self.mask_dir,'atlases','HarvardOxford','HarvardOxford-cort-maxprob-thr0-2mm.nii.gz')
        self.ho_subcortical = os.path.join(self.mask_dir,'atlases','HarvardOxford','HarvardOxford-sub-maxprob-thr0-2mm.nii.gz')
        
        # load cortical
        ho_cortical = np.array(nib.load(self.ho_cortical).get_data(),dtype=float)
        
        # load subcortical
        ho_subcortical = np.array(nib.load(self.ho_cortical).get_data(),dtype=float)
        
        shell()
        # where cortical == 0, make the label subcortical
        
        ho = np.where(ho_cortical==0,ho_subcortical,ho_cortical)
        

    def roy_rsa_letters(self,task='rsa'):
        # Use output from the rsa_letters analysis
        # For each letter, extract the t-stats from all voxels
        
        mask = np.array(nib.load(self.mni).get_data(),dtype=bool) # boolean mask
        labels = np.array(nib.load(self.mni_labels).get_data(),dtype=float) # labels

        letters = [
            'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
            'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'
        ] 
        
        DFOUT = pd.DataFrame() # output dataframe for all subjects
        for s,subject in enumerate(self.subjects):
            for sess,session in enumerate(self.sessions):
                # path to first level feat directory
                this_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_letters_sub-{}_{}.feat'.format(task,subject,session),'stats')
                
                # 52 EVs alphabet in black, then alphabet in color
                for cope in np.arange(1,53):
                    this_df = pd.DataFrame() # temporary DF for concatenation
                    BOLD = os.path.join(this_path,'tstat{}.nii.gz'.format(cope)) 

                    nii = nib.load(BOLD).get_data()[mask] # flattens
                    # columns
                    this_df['subject']  = np.repeat(subject,len(nii))
                    this_df['session']  = np.repeat(sess+1,len(nii))
                    this_df['ev']       = np.repeat(cope,len(nii))
                    this_df['tstat']    = nii
                    this_df['brain_labels']  = labels[mask] # flatten
                    # make letter and color_condition columns based on EV numbers
                    if cope < 27:
                        color_condition = 'black' # first are all black
                    else:
                        color_condition = 'color'    
                    this_df['letter'] = np.repeat(letters[cope-1],len(nii))
                    this_df['color_condition'] = np.repeat(color_condition, len(nii))
                    print('cope{} letter={} {}'.format(cope, letters[cope-1], color_condition))
                    
                    DFOUT = pd.concat([DFOUT,this_df],axis=0)
        os.path.join(self.higher_level_dir,'task-{}'.format(task),'roy_task-{}_letters.csv'.format(task))
        # DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'roy_task-{}_letters.csv'.format(task)),sep='/t')
        
        
        
        print('success: roy_rsa_letters')
        
        
        
        