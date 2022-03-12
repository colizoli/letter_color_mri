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
        self.ho_cortical_labels = os.path.join(self.mask_dir,'atlases','HarvardOxford','HarvardOxford-{}-maxprob-thr0-2mm.nii.gz'.format('cort'))
        self.ho_cortical_prob = os.path.join(self.mask_dir,'atlases','HarvardOxford','HarvardOxford-{}-prob-2mm.nii.gz'.format('cort'))
        
        # write unix commands to job to run in parallel
        self.higher_level_job_path = os.path.join(self.analysis_dir,'jobs','job_higher_level.txt')
        if not os.path.exists(self.higher_level_job_path):
            self.higher_level_job_path = open(self.higher_level_job_path, "w")
            self.higher_level_job_path.write("#!/bin/bash\n")
            self.higher_level_job_path.close()
        
        self.letters = [
            'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
        ]
       
    def rsa_letters_combine_colors(self,):      
    # concatenates all subjects colors files into a single one
        
        DFOUT = pd.DataFrame()
        for s,subject in enumerate(self.subjects):
            colors = pd.read_csv(os.path.join(self.deriv_dir,'sub-{}'.format(subject),'sub-{}_colors.tsv'.format(subject)),sep='\t')
            # concate
            DFOUT = pd.concat([DFOUT,colors],axis=0)
        
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format('rsa'),'task-rsa_trained_colors.tsv'), sep='\t')
        print('success: rsa_letters_combine_colors')
    
    
    
    def labels_harvard_oxford(self,):
        # outputs the labels and probabilities as a vector
        # harvard oxford cortical and subcortical separately (note subcortical probabilities are only on 1 hemisphere)
        
        for atlas in ['cort',]:
            DFOUT = pd.DataFrame() 
            # load 
            mask        = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=bool) # binary mask
            labels      = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=float)
            probability = np.array(nib.load(self.ho_cortical_prob).get_data(),dtype=float) # ROIs in time axis
            
            # load subcortical
            DFOUT['{}_labels'.format(atlas)]  = labels[mask] # flatten
            
            # loop through all labels and extract probabilities
            for this_label in np.arange(np.max(labels)):                
                roi = probability[:,:,:,int(this_label)] # time axis for current ROI
                DFOUT['{}_{}'.format(atlas,int(this_label+1))]  = roi[mask] # save probabilities as new column
                        
            DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format('rsa'),'harvard_oxford_{}_probabilities_flattened.csv'.format(atlas)),sep='\t')
        print('success: labels_harvard_oxford')
        
            
    def rsa_letters_ev_conditions(self,task='rsa'):
        # outputs a dataframe with the letter and color conditions (general) for the EVs in the rsa-letters analysis
        # a-z black, then a-z color, 53rd EV was oddballs (nuisance regressor)

        #############################
        # output ev-number to letter-color condition file
        # make letter and color_condition columns based on EV numbers
        LCDF = pd.DataFrame()
        LCDF['ev'] = np.arange(1,53)
        LCDF['letter'] = self.letters*2        
        LCDF['color_condition'] = np.concatenate([np.repeat('black',len(letters)),np.repeat('color',len(letters))])
        LCDF.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_letters_ev_conditions.csv'.format(task)),sep='\t')
        print('success: rsa_letters_ev_conditions')
        
        
    def roy_rsa_letters(self,task='rsa'):
        # Use output from the rsa_letters analysis
        # For each letter, extract the t-stats from all voxels
        
        # using Harvard Oxford cortical
        mask = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=bool) # boolean mask
        labels = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=float) # labels
        
        DFOUT = pd.DataFrame() # output tstat dataframe for all subjects
        for s,subject in enumerate(self.subjects):
            for sess,session in enumerate(self.sessions):
                # path to first level feat directory
                this_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_letters_sub-{}_{}.feat'.format(task,subject,session),'stats')
                
                # 52 EVs alphabet in black, then alphabet in color
                for cope in np.arange(1,53):
                    this_df = pd.DataFrame() # temporary DF for concatenation
                    BOLD = os.path.join(this_path,'tstat{}.nii.gz'.format(cope)) 
                    # statistic
                    nii = nib.load(BOLD).get_data()[mask] # flattens
                    # columns
                    this_df['subject']  = np.repeat(subject,len(nii))
                    this_df['session']  = np.repeat(sess+1,len(nii))
                    this_df['ev']       = np.repeat(cope,len(nii))
                    this_df['tstat']    = nii
                    this_df['brain_labels']  = labels[mask] # flatten

                    # concat data frames
                    DFOUT = pd.concat([DFOUT,this_df],axis=0)
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'roy_task-{}_letters.csv'.format(task)),sep='\t')
        print('success: roy_rsa_letters')
        
    def kelly_rsa_letters(self,task='rsa'):
        # Use output from the rsa_letters analysis
        # For each letter, extract the z-stats from OCCIPITAL LOBE
        # use occipital lobe mask from MNI atlas = label #5
        
        N = 5 # label number
        mni_labels = np.array(nib.load(self.mni_labels).get_data(),dtype=float)
        mask = np.array(np.where(mni_labels==N,1,0),dtype=bool) # binary mask, this roi
                
        DFOUT = pd.DataFrame() # output tstat dataframe for all subjects
        for s,subject in enumerate(self.subjects):
            for sess,session in enumerate(self.sessions):
                # path to first level feat directory
                this_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_letters_sub-{}_{}.feat'.format(task,subject,session),'stats')
                
                # 52 EVs alphabet in black, then alphabet in color
                for cope in np.arange(1,53):
                    this_df = pd.DataFrame() # temporary DF for concatenation
                    BOLD = os.path.join(this_path,'zstat{}.nii.gz'.format(cope)) 
                    # statistic
                    nii = nib.load(BOLD).get_data()[mask] # flattens
                    # columns
                    this_df['subject']  = np.repeat(subject,len(nii))
                    this_df['session']  = np.repeat(sess+1,len(nii))
                    this_df['ev']       = np.repeat(cope,len(nii))
                    this_df['zstat']    = nii
                    this_df['brain_labels']  = np.repeat(N,len(nii)) # flatten

                    # concat data frames
                    DFOUT = pd.concat([DFOUT,this_df],axis=0)
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'kelly_task-{}_letters.csv'.format(task)),sep='\t')
        print('success: kelly_rsa_letters')    
        
        