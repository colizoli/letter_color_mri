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
            self.higher_level_job = open(self.higher_level_job_path, "w")
            self.higher_level_job.write("#!/bin/bash\n")
            self.higher_level_job.close()
        
        self.letters = [
            'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
        ]
       
    def localizers_randomise_input(self,task):
        # Localizers (letters, colors), cope1 was contrast of interest
        # For randomise input, we need the 3D cope1 images stacked for all subjects in the 4th dimension
        # i.e., MNI x #subjects
        
        if task == 'letters':
            mask = os.path.join(self.mask_dir,'occipital_fusiform_temporal_occipital_fusiform_cortex_LH.nii.gz') # search only within these voxels
        elif task == 'colors':
            mask =  os.path.join(self.mask_dir,'occipital_temporal_occipital_fusiform_cortex.nii.gz')
        
        # output path MNI x #subjects
        outFile = os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_cope1.nii.gz'.format(task))
        # Load MNI brain and save headers
        mni = nib.load(self.mni)
        # initialize empty 4D array
        COPE = np.zeros((mni.get_data().shape[0],mni.get_data().shape[1],mni.get_data().shape[2],len(self.subjects)))
        for s,subject in enumerate(self.subjects):
            C1_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}.feat'.format(task,subject),'stats','cope1.nii.gz')
            C1 = nib.load(C1_path).get_data()
            COPE[:,:,:,s] = C1
            print(subject)
        # save as nifti file with correct header and data type        
        outData = nib.Nifti1Image(COPE, affine=mni.affine, header=mni.header) # pass affine and header from last MNI image
        outData.set_data_dtype(np.float32)
        nib.save(outData, outFile)
        
        # open higher level job and write command as new line
        randomise_output = os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_cope1'.format(task))
        cmd = 'randomise -i {} -o {} -m {} -1 -x -c 3.1'.format(outFile,randomise_output,mask) #outFile is from concatenation
        self.higher_level_job = open(self.higher_level_job_path, "a") # append is important, not write
        self.higher_level_job.write(cmd)   # feat command
        self.higher_level_job.write("\n\n")  # new line
        self.higher_level_job.close()
        print('success: localizers_randomise_input')
        
    def rsa_letters_ev_conditions(self,task='rsa'):
        # outputs a dataframe with the letter and color conditions (general) for the EVs in the rsa-letters analysis
        # a-z black, then a-z color, 53rd EV was oddballs (nuisance regressor)

        #############################
        # output ev-number to letter-color condition file
        # make letter and color_condition columns based on EV numbers
        LCDF = pd.DataFrame()
        LCDF['ev'] = np.arange(1,53)
        LCDF['letter'] = self.letters*2        
        LCDF['trial_type_color'] = np.concatenate([np.repeat('black',len(self.letters)),np.repeat('color',len(self.letters))])
        LCDF.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_letters_ev_conditions.tsv'.format(task)),sep='\t')
        print('success: rsa_letters_ev_conditions')
        
    def rsa_letters_conditions(self,task='rsa'):      
        # concatenates all subjects colors files into a single one
        
        # load evs conditions file
        LCDF = pd.read_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_letters_ev_conditions.csv'.format(task)),sep='\t')
        try:
            LCDF = LCDF.loc[:, ~LCDF.columns.str.contains('^Unnamed')]
        except:
            pass
        
        DFOUT = pd.DataFrame()
        for s,subject in enumerate(self.subjects):
            events = pd.read_csv(os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}_{}_events.tsv'.format(task,subject,'ses-01')),sep='\t')
            try:
                events = events.loc[:, ~events.columns.str.contains('^Unnamed')]
            except:
                pass
            # drop trial-wise information
            events = events[events['oddball']==0]
            colors = events.get(['subject','letter', 'trial_type_letter', 'trial_type_color','r','g', 'b',])
            # drop duplicates
            colors = colors.drop_duplicates()
            # a-z black, then a-z color
            colors = colors.sort_values(by=['trial_type_color','letter']).reset_index()
            # merge with ev condition file
            colors = pd.merge(colors, LCDF, on=["letter", "trial_type_color"])
            # concate
            DFOUT = pd.concat([DFOUT,colors],axis=0)
        
        DFOUT = DFOUT.drop(['index']) # drop the index col
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format('rsa'),'task-rsa_letters_conditions.tsv'), sep='\t')
        print('success: rsa_letters_combine_colors')
    
    def labels_harvard_oxford(self,):
        # outputs the labels and probabilities as a vector
        # harvard oxford (sub)cortical atlas
        # note that subcortical atlas has left/right hemispheres as separate ROIs
        
        for atlas in ['cort',]:
            DFOUT = pd.DataFrame() 
            # load 
            mask        = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=bool) # binary mask
            labels      = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=float) # labels as integers in mask
            probability = np.array(nib.load(self.ho_cortical_prob).get_data(),dtype=float) # ROIs in time axis
            
            # add column labels
            DFOUT['{}_labels'.format(atlas)]  = labels[mask] # flatten
            
            # loop through all labels and extract probabilities
            for this_label in np.arange(np.max(labels)):                
                roi = probability[:,:,:,int(this_label)] # time axis for current ROI
                DFOUT['{}_{}'.format(atlas,int(this_label+1))]  = roi[mask] # save probabilities as new column
            DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format('rsa'),'harvard_oxford_{}_probabilities_flattened.tsv'.format(atlas)),sep='\t')
        print('success: labels_harvard_oxford')
    
    def probabilities_emotion_rois(self,):
        # outputs the probabilities of each brain region as columns
        # used a combination of the harvard oxford cortical and subcortical atlases
        # there are many overlapping voxels, i.e., voxels that have non-zero probabilities in 2+ ROIs
        
        roi_path = os.path.join(self.mask_dir,'hilde_brain_regions')
        DFOUT = pd.DataFrame() 
        for roi in os.listdir(roi_path):
            # load             
            mask        = np.array(nib.load(os.path.join(self.mask_dir,'emotion_brain_regions_mask.nii.gz')).get_data(),dtype=bool) # binary mask
            probability = np.array(nib.load(os.path.join(roi_path,roi)).get_data(),dtype=float) # ROIs in time axis
            # save probabilities as new column
            DFOUT[roi[:-7]]  = probability[mask] # (remove trailing '.nii.gz')
            DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format('rsa'),'hilde_brain_region_probabilities.tsv'),sep='\t')
            print(roi)
        print('success: labels_harvard_oxford')
        
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
                    this_df['mask_idx']      = np.arange(len(nii))

                    # concat data frames
                    DFOUT = pd.concat([DFOUT,this_df],axis=0)
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'roy_task-{}_letters.tsv'.format(task)),sep='\t')
        print('success: roy_rsa_letters')
    
    
    def hilde_rsa_letters(self,task='rsa'):
        # Use output from the rsa_letters analysis
        # For each letter, extract the t-stats from all voxels in "emotional ROIs"
        
        # using Harvard Oxford cortical
        emotion_rois = os.path.join(self.mask_dir,'emotion_brain_regions_mask.nii.gz')
        mask = np.array(nib.load(emotion_rois).get_data(),dtype=bool) # boolean mask
        # labels = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=float) # labels
        
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
                    # this_df['brain_labels']  = labels[mask] # flatten
                    this_df['mask_idx']      = np.arange(len(nii))

                    # concat data frames
                    DFOUT = pd.concat([DFOUT,this_df],axis=0)
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'hilde_task-{}_letters.tsv'.format(task)),sep='\t')
        shell()
        print('success: hilde_rsa_letters')
        
        
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
                    this_df['mask_idx']      = np.arange(len(nii))
                    
                    # concat data frames
                    DFOUT = pd.concat([DFOUT,this_df],axis=0)
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'kelly_task-{}_letters.tsv'.format(task)),sep='\t')
        print('success: kelly_rsa_letters')    
        
    def timeseries_trials_rsa(self,kernel,task='rsa'):
        # For each trial in RSA task, extract time series data
        # Timeseries with length = kernel (#samples)
        # Input is the same input to the first level analysis rsa_letters
        
        # using Harvard Oxford cortical
        mask = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=bool) # boolean mask
        labels = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=float) # labels
        
        DFOUT = pd.DataFrame() # output tstat dataframe for all subjects
        for s,subject in enumerate(self.subjects):
            for sess,session in enumerate(self.sessions):
                # need trial-wise information with event onsets
                events = pd.read_csv(os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}_{}_events.tsv'.format(task,subject,session)),sep='\t')
                events = events.loc[:, ~events.columns.str.contains('^Unnamed')] # drop unnamed columns
                BOLD_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}_{}_bold_mni.nii.gz'.format(task,subject,session)) 
                BOLD = nib.load(BOLD_path).get_data() # numerical data 4D 
                
                # convert stimulus onset to TR, round to nearest TR
                events['nearest_TR'] = np.round(events['onset']/float(self.TR))
                TRs = np.array(events['nearest_TR'])
                
                this_df = [] # temporary list for concatenation
                # loop through trials in events dataframe
                for trial in np.arange(events.shape[0]):
                    start_tr = int(TRs[trial])-1 # start at TR-1 for indexing 
                    stop_tr = int(TRs[trial])-1+kernel # extract time series with length = kernel
                    this_series = BOLD[:,:,:,start_tr:stop_tr] # minus 1 for index of TR   
                    this_series_flat =  this_series[mask]  # flatten 3D brain to 1D
                    # trials x voxels x kernel
                    this_df.append(this_series_flat)
                
                this_df = np.array(this_df) # trials x voxels x kernel
                # output as numpy array
                
                # grouped = events.groupby(['letter','trial_type_color'])
                shell()
        print('success: timeseries_trials_rsa')
    
    def timeseries_letters_rsa(self,kernel,task='rsa'):
        # For each letter-color condition in RSA task, extract time series data (a-z black, then a-z color)
        # Timeseries with length = kernel (#samples)
        # Input is the same input to the first level analysis rsa_letters
        
        # using Harvard Oxford cortical
        mask = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=bool) # boolean mask
        labels = np.array(nib.load(self.ho_cortical_labels).get_data(),dtype=float) # labels
        
        DFOUT = pd.DataFrame() # output tstat dataframe for all subjects
        for s,subject in enumerate(self.subjects):
            for sess,session in enumerate(self.sessions):
                # need trial-wise information with event onsets
                events = pd.read_csv(os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}_{}_events.tsv'.format(task,subject,session)),sep='\t')
                events = events.loc[:, ~events.columns.str.contains('^Unnamed')] # drop unnamed columns
                BOLD_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}_{}_bold_mni.nii.gz'.format(task,subject,session)) 
                BOLD = nib.load(BOLD_path).get_data() # numerical data 4D 
                
                # convert stimulus onset to TR, round to nearest TR
                events['nearest_TR'] = np.round(events['onset']/float(self.TR))
                
                # loop through trials in events dataframe
                for C in np.unique(events['trial_type_color']):
                    # loop letters
                    for L in self.letters:
                        # check if letter-color condition exists, get only those trials in events file
                        this_ev = events[(events['letter']==L) & (events['trial_type_color']==C)]
                        TRs = np.array(this_ev['nearest_TR'])
                        
                        if not this_ev.empty:
                            # CUT OUT TIME SERIES 
                            this_ts_df = [] # temporary DF for letter-color condition all trials
                            # go through onsets only in current letter-color condition
                            for trial in np.arange(this_ev.shape[0]):
                                start_tr = int(TRs[trial])-1 # start at TR-1 for indexing 
                                stop_tr = int(TRs[trial])-1+kernel # extract time series with length = kernel
                                this_series = BOLD[:,:,:,start_tr:stop_tr] # minus 1 for index of TR   
                                this_series_flat =  this_series[mask]  # flatten 3D brain to 1D
                                this_ts_df.append(this_series_flat)
                            # take mean timeseries kernel of all trials in current letter-color condition
                            nii = np.mean(np.array(this_ts_df),axis=0)
                            this_df = pd.DataFrame() # temporary DF for concatenation, this subject, session
                            # trials x voxels x kernel
                            this_df['subject']  = np.repeat(subject,len(nii))
                            this_df['session']          = np.repeat(sess+1,len(nii))
                            this_df['letter']           = np.repeat(L,len(nii))
                            this_df['trial_type_color'] = np.repeat(C,len(nii))
                            this_df['brain_labels']     = labels[mask] # flatten
                            this_df['mask_idx']         = np.arange(len(nii))
                        
                            nii = pd.DataFrame(nii)
                            # concat data frames
                            this_df = pd.concat([nii,this_df],axis=1) # add kernel as columns in front 
                            DFOUT = pd.concat([this_df,DFOUT],axis=0) # add entire DF as rows to bottom
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_letters_timeseries.tsv'.format(task)),sep='\t')
        print('success: timeseries_letters_rsa')