#!/usr/bin/env python
# encoding: utf-8

"""
letter_color_housekeeping.py
Created by O.Colizoli, June-2025
Last update: 09-01-2025
Python version 3.9

Change filenames, move files, etc.
"""

import os, subprocess, sys, glob
import shutil as sh
import pandas as pd
import numpy as np
import json
from datetime import datetime
from IPython import embed as shell # for Oly's debugging only

class housekeeping_class(object):
    """
    Define a class for the housekeeping of the MRI data set.
    
    Args:
        subject (str): Subject number main experiment.
        bids_dir (str): Path to the DICOM files for the MRI acquisition.
        
    Attributes:
        subject (str): Subject number main experiment.
        bids_dir (str): Path to the NIFTI files after converstion from DICOM.
    
    Methods:
        __init__(subject, subject, bids_dir)
        rename_behav_logfiles()
    """
    
    
    def __init__(self, subject, bids_dir, deriv_dir):        
        self.subject        = 'sub-'+str(subject)
        self.bids_dir       = str(bids_dir)
        self.deriv_dir       = str(deriv_dir)
        # path to fmriprep output (preprocessed fmri data)
        self.fmriprep_dir = os.path.join(self.deriv_dir, 'fmriprep')
    
    
    def rename_behav_logfiles(self):
        """Rename the behavioural log files to be bids compliant: e.g., sess-2 -> ses-02.
        """
        old_ses = ['sess-1','sess-2']
        
        for s,new_ses in enumerate(['ses-mri01', 'ses-mri02']):
            ###################
            # RENAME 'behav' 'sess-1' to 'ses-01' in sub-xxx/ses-mri01/behav
            ###################
            dir_path = os.path.join(self.bids_dir, self.subject, new_ses, 'behav')
            
            for f in os.listdir(dir_path):
                if old_ses[s] in f:
                    f_new = f.replace(old_ses[s], new_ses)
                    os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                    print('old={} , new={}'.format(f,f_new))
                                        
        print('success: rename_behav_logfiles')
    

        
    #
    # def loc_rename_bold(self):
    #     """Match the LOC1 and LOC2 bold acquisition (nifti) files to the letters and colors localizer events.
    #         Rename nifti files to include localizer name.
    #     """
    #     T = 6 # number of trailing timestamp characters in events file names
    #
    #     for session in ['ses-mri01', 'ses-mri02']:
    #         for preprocessed_tag in ['space-MNI152NLin6Asym_desc-preproc_bold', 'space-T1w_desc-preproc_bold']:
    #             ###################
    #             # check acquisition time in JSON files
    #             ###################
    #             # JSON paths
    #             loc1_path = os.path.join(self.bids_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC01_dir-AP_bold.json'.format(self.subject, session))
    #             loc2_path = os.path.join(self.bids_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC02_dir-AP_bold.json'.format(self.subject, session))
    #
    #             # LOC1 nifti
    #             with open(loc1_path) as f:
    #                 loc1_json = json.load(f)
    #                 loc1_acq = loc1_json['AcquisitionTime']
    #                 loc1_acq = os.path.splitext(loc1_acq)[0] # drop extension
    #
    #             # LOC2 nifti
    #             with open(loc2_path) as f:
    #                 loc2_json = json.load(f)
    #                 loc2_acq = loc2_json['AcquisitionTime']
    #                 loc2_acq = os.path.splitext(loc2_acq)[0] # drop extension
    #
    #             # letters events
    #             for f1 in os.listdir(os.path.join(self.bids_dir, self.subject, session, 'func')):
    #                 if ('task-letters' in f1) and ('events' in f1):
    #                     fname, extension = os.path.splitext(f1) # split name from extension
    #                     letters_acq = fname[-T:]
    #
    #             # colors events
    #             for f2 in os.listdir(os.path.join(self.bids_dir, self.subject, session, 'func')):
    #                 if ('task-colors' in f2) and ('events' in f2):
    #                     fname, extension = os.path.splitext(f2) # split name from extension
    #                     colors_acq = fname[-T:]
    #
    #             ### MATCHING ###
    #             letters_first = int(letters_acq) < int(colors_acq) # if letters are first, acquisition time is less
    #             loc1_first = datetime.strptime(loc1_acq, '%H:%M:%S') < datetime.strptime(loc2_acq, '%H:%M:%S') # if loc1 is first, acquisition time is less
    #
    #             if letters_first:
    #                 if loc1_first:
    #                     # letters == loc1
    #                     # colors = loc2
    #                     # NIFTI paths
    #                     old = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC01_dir-AP_{}.nii.gz'.format(self.subject, session, preprocessed_tag))
    #                     new = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC01letters_dir-AP_{}.nii.gz'.format(self.subject, session, preprocessed_tag))
    #                     os.rename(old, new) # old,new
    #                     print('old={} , new={}'.format(old, new))
    #                     old = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC02_dir-AP_{}.nii.gz'.format(self.subject, session, preprocessed_tag))
    #                     new = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC02colors_dir-AP_{}.nii.gz'.format(self.subject, session, preprocessed_tag))
    #                     os.rename(old, new) # old,new
    #                     print('old={} , new={}'.format(old, new))
    #                 else:
    #                     # colors = loc1
    #                     # letters  = loc2
    #                     old = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC01_dir-AP_{}.nii.gz'.format(self.subject, session, preprocessed_tag))
    #                     new = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC01colors_dir-AP_{}.nii.gz'.format(self.subject, session, preprocessed_tag))
    #                     os.rename(old, new) # old,new
    #                     print('old={} , new={}'.format(colors_old, colors_new))
    #                     old = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC02_dir-AP_{}.nii.gz'.format(self.subject, session, preprocessed_tag))
    #                     new = os.path.join(self.fmriprep_dir, self.subject, session, 'func', '{}_{}_task-cmrr2isomb4TR1500LOC02letters_dir-AP_{}.nii.gz'.format(self.subject, session, preprocessed_tag))
    #                     os.rename(old, new) # old,new
    #                     print('old={} , new={}'.format(old, new))
    #
    #             print('{} {} loc1_first={} letters_first={}'.format(self.subject, session, letters_first, loc1_first))
    #     print('success: loc_rename_bold')
    #
        
            
    # def remove_timestamps_logfiles(self, ):
    #     """Remove the trailing timestamps from the psychopy output: last 16 characters _00200312-125726.
    #
    #     Args:
    #         behav_sess (str): The session string within the behavioral logfile names.
    #     """
    #     T = 16 # length of trailing characters to remove
    #
    #     ###################
    #     # behav folder
    #     ###################
    #     dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'behav')
    #     # loops through functional files
    #     for f in os.listdir(dir_path):
    #         if ('task-iconic' in f) or ('prefs' in f) or ('task-consistency' in f):
    #             Fname, extension = os.path.splitext(f) # split name from extension
    #             F_new = Fname[:-T]+extension # remove T trailing T characters, add extension back
    #             os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
    #             print('old={} , new={}'.format(f,f_new))
    #
    #     ###################
    #     # func folder
    #     ###################
    #     dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
    #     # loops through functional files
    #     for f in os.listdir(dir_path):
    #         fname, extension = os.path.splitext(f) # split name from extension
    #         f_new = fname[:-T]+extension # remove T trailing T characters, add extension back
    #         os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
    #         print('old={} , new={}'.format(f,f_new))
    #
    #     print('success: remove_timestamps_logfiles')
    #
    # def copy_logfiles(self, behav_sess):
    #     """Copy events and log files into the bids_raw directory with the NIFTI files.
    #
    #     Args:
    #         behav_sess (str): The session string within the behavioral logfile names.
    #     """
    #     ###################
    #     # copy colors file
    #     ###################
    #     src = os.path.join(self.raw_dir, 'logfiles', self.subject, '{}_colors.tsv'.format(self.subject))
    #     dst = os.path.join(self.raw_dir, 'nifti', self.subject, '{}_colors.tsv'.format(self.subject))
    #     sh.copyfile(src, dst)
    #     print('colors: src={} , dst={}'.format(src,dst))
    #
    #     ###################
    #     # copy 'behav' folder
    #     ###################
    #     src = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'behav')
    #     dst = os.path.join(self.raw_dir, 'nifti', self.subject, self.session, 'behav')
    #     sh.copytree(src, dst)
    #     print('behav: src={} , dst={}'.format(src,dst))
    #
    #     ###################
    #     # copy mri event files into existing 'func' folder
    #     ###################
    #     # loop through functional files
    #     for f in os.listdir(os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')):
    #         src = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func', f)
    #         dst = os.path.join(self.raw_dir, 'nifti', self.subject, self.session, 'func', f)
    #         sh.copy(src, dst)
    #     print('func: src={} , dst={}'.format(src,dst))
    #
    #     print('success: copy_logfiles')

        

