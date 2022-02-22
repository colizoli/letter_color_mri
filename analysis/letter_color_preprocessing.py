#!/usr/bin/env python
# encoding: utf-8

"""
letter_color_preprocessing.py
Created by O.Colizoli, June-2019
Last update: 11-06-2019
Python version 3.6

The following packages need to be installed because they are called from the command line, but not imported:
dcm2niix (conda install -c conda-forge dcm2niix)
dcm2bids (python -m pip install dcm2bids)
fsl

"""

import os, subprocess, sys
import shutil as sh
import nibabel as nib
import pandas as pd
from IPython import embed as shell # for Oly's debugging only

class preprocess_class(object):
    def __init__(self, subject, mri_subject, session, analysis_dir, source_dir, raw_dir, deriv_dir, mask_dir, template_dir, timing_files_dir):        
        self.subject        = 'sub-'+str(subject)
        self.mri_subject    = 'sub-'+str(mri_subject)
        self.session        = str(session)
        self.analysis_dir   = str(analysis_dir)
        self.source_dir     = str(source_dir)   # dicoms
        self.raw_dir        = str(raw_dir)      # raw niftis
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.timing_files_dir = str(timing_files_dir)
            
        if not os.path.isdir(self.raw_dir):
            os.mkdir(self.raw_dir)
        if not os.path.isdir(self.deriv_dir):
            os.mkdir(self.deriv_dir)
        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
        if not os.path.isdir(self.timing_files_dir):
            os.mkdir(self.timing_files_dir)    

    
    def dicom2bids(self, ):
        """ Converts MRI data from dicom to nifti format, then structures them according to the bids standard. 
        Takes about 15 minutes per session.
        Requirements:
            dcm2niix (conda install -c conda-forge dcm2niix)
            dcm2bids (python -m pip install dcm2bids) 
        """
        # for coverting from source directories to bids file names
        if self.session == 'ses-01':
            mri_session =  'ses-mri01'
        else:
            mri_session = 'ses-mri02'
        
        config = 'convert2bids_letter-color_main_{}.json'.format(self.subject)
        
        DICOM_DIR = os.path.join(self.source_dir, self.mri_subject, mri_session) #/Users/olympiacolizoli/Aeneas/mountpoint6/raw
        CONFIG_FILE  = os.path.join(self.analysis_dir, 'config', config)
        OUTPUT_DIR = os.path.join(self.raw_dir, 'nifti')

        # Create folder structure if not yet existing
        if not os.path.isdir(os.path.join(self.raw_dir, 'nifti')):
            os.mkdir(os.path.join(self.raw_dir, 'nifti'))

        cmd = 'dcm2bids -d {} -p {} -s {} -c {} -o {}'.format(DICOM_DIR, self.subject, self.session, CONFIG_FILE, OUTPUT_DIR)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        print('converting sub-{} convert2bids_letter-color_main_{}.json'.format(self.mri_subject,self.subject))
        print('success: dicom2bids')
    
    def rename_logfiles(self, behav_sess):
        """ Renames the log files to match MRI files: e.g., sess-2 -> ses-02
        """
        ###################
        # RENAME 'behav' 'sess-1' to 'ses-01' in logfiles/behav folder before copying
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'behav')
        for F in os.listdir(dir_path):
            if behav_sess in F:
                F_new = F.replace(behav_sess, self.session)
                os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                print('old={} , new={}'.format(F,F_new))
        
        
        ###################
        # RENAME 'func' 'sess-1' to 'ses-01' in logfiles/func folder before copying
        # then runs 'run-1' to 'run-01' etc
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
        behav_runs = ['run-1','run-2','run-3','run-4']
        bids_runs = ['run-01','run-02','run-03','run-04']
        # loops through functional files
        for F in os.listdir(dir_path):
            # replace sessions
            if behav_sess in F:
                F_new = F.replace(behav_sess, self.session)
                os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                print('old={} , new={}'.format(F,F_new))
        # restart loop because new session changed names
        for F in os.listdir(dir_path):
            # replaces runs
            for r,this_run in enumerate(behav_runs):
                if this_run in F:
                    F_new = F.replace(this_run, bids_runs[r])
                    os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                    print('old={} , new={}'.format(F,F_new))
        
        ###################
        # RENAME localizer runs only: remove '_run-01' because only 1 run of each
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
        # loops through functional files
        for F in os.listdir(dir_path):
            # remove run from localizer event file names
            if ('task-colorsloc' in F) or ('task-lettersloc' in F):
                F_new = F.replace('_run-01', '')
                os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                print('old={} , new={}'.format(F,F_new))
                try: # some files have runs 5 or 6
                    F_new = F.replace('_run-5', '')
                    os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                    F_new = F.replace('_run-6', '')
                    os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                except:
                    pass
                    
        print('success: rename_logfiles')
    
    def remove_timestamps_logfiles(self, behav_sess):
        """ removes the trailing timestamps from the psychopy output
         e.g. last 16 characters _00200312-125726
        """
        T = 16 # length of trailing characters to remove
        
        ###################
        # behav folder
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'behav')
        # loops through functional files
        for F in os.listdir(dir_path):
            if ('task-iconic' in F) or ('prefs' in F) or ('task-consistency' in F):
                Fname, extension = os.path.splitext(F) # split name from extension
                F_new = Fname[:-T]+extension # remove T trailing T characters, add extension back
                os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                print('old={} , new={}'.format(F,F_new))
                
        ###################
        # func folder
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
        # loops through functional files
        for F in os.listdir(dir_path):
            Fname, extension = os.path.splitext(F) # split name from extension
            F_new = Fname[:-T]+extension # remove T trailing T characters, add extension back
            os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
            print('old={} , new={}'.format(F,F_new))
        
        print('success: remove_timestamps_logfiles')
        
    def copy_logfiles(self, behav_sess):
        """ copies events and log files into the bids_raw directory with the nifti files
        """
        ###################
        # copy colors file
        ###################
        src = os.path.join(self.raw_dir, 'logfiles', self.subject, '{}_colors.tsv'.format(self.subject))
        dst = os.path.join(self.raw_dir, 'nifti', self.subject, '{}_colors.tsv'.format(self.subject))
        sh.copyfile(src, dst)
        print('colors: src={} , dst={}'.format(src,dst))

        ###################
        # copy 'behav' folder
        ###################
        src = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'behav')
        dst = os.path.join(self.raw_dir, 'nifti', self.subject, self.session, 'behav')
        sh.copytree(src, dst)
        print('behav: src={} , dst={}'.format(src,dst))
        
        ###################
        # copy mri event files into existing 'func' folder
        ###################
        # loop through functional files
        for F in os.listdir(os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')):
            src = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func', F)
            dst = os.path.join(self.raw_dir, 'nifti', self.subject, self.session, 'func', F)
            sh.copy(src, dst)
        print('func: src={} , dst={}'.format(src,dst))
        
        print('success: copy_logfiles')
        
    def housekeeping(self, ):
        """ Renames then copies events and log files into the bids_raw directory with the nifti files
        """
        if self.session == 'ses-01':
            behav_sess = 'sess-1'
        else:
            behav_sess = 'sess-2'
        
        self.rename_logfiles(behav_sess)    # rename session and runs
        ## WARNING!! don't remove timestamps twice on same participant
        self.remove_timestamps_logfiles(behav_sess)  # removes the trailing timestamps from the psychopy output
        self.copy_logfiles(behav_sess)      # copy logfiles to bids_raw folder with mri data
        print('success: housekeeping')
    
    def raw_copy(self,):
        """ copy from bids_raw directory into derivaties folder for further processing/analysis.
        All further analysis should be run within the derivatives folder.
        """
        if os.path.isdir(os.path.join(self.deriv_dir, self.subject, self.session)):
            print('Target folder already exists. No files will be copied. Delete existing session from derivatives if intended.')
            return
        else:
            # Copy nifti
            src = os.path.join(self.raw_dir, 'nifti', self.subject, self.session)
            dst = os.path.join(self.deriv_dir, self.subject, self.session)
            sh.copytree(src, dst)
            # copy colors file
            src = os.path.join(self.raw_dir, 'nifti', self.subject, '{}_colors.tsv'.format(self.subject))
            dst = os.path.join(self.deriv_dir, self.subject, '{}_colors.tsv'.format(self.subject))
            sh.copyfile(src, dst)
            print('colors: src={} , dst={}'.format(src,dst))
        print('success: raw_copy')

    def bet_brains_T1(self, postfix='brain'):
        """ Runs brain extraction on all T1s.
        Always check visually in fsleyes.
        """
    
        inFile = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w.nii.gz'.format(self.subject, self.session))
        outFile = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w_{}.nii.gz'.format(self.subject, self.session, postfix))
        # bet inFile outFile
        cmd = 'bet {} {}'.format(inFile, outFile)
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)

        # reorient to mni space if necessary (sometimes turned around)
        # cmd = 'fslreorient2std {} {}'.format(outFile,outFile)
        # proc = subprocess.call( cmd, shell=True, bufsize=0,) # COMMAND LINE
        print('success: bet_brains_T1')
    
    def bet_brains_fmap(self, postfix='brain'):
        # Runs brain extraction on magnitude images of field maps
        # This needs to be 'tight', meaning it is better to exclude brain voxels than include noisy non-brain voxels
        # Important! Always check visually in fsleyes

        for run in ['01','02','03']: # first 2 are magnitude, 3rd is phase
            inFile = os.path.join(self.deriv_dir,  self.subject,  self.session, 'fmap', '{}_{}_run-{}_fmap.nii.gz'.format(self.subject,self.session,run))
            outFile = os.path.join(self.deriv_dir,  self.subject,  self.session, 'fmap', '{}_{}_run-{}_fmap_{}.nii.gz'.format(self.subject,self.session,run,postfix))
            # bet inFile outFile
            cmd = 'bet {} {} -f 0.6'.format(inFile, outFile) 
            print(cmd)
            results = subprocess.call(cmd, shell=True, bufsize=0)
    
            # reorient to mni space if necessary (sometimes turned around)
            # cmd = 'fslreorient2std {} {}'.format(outFile,outFile)
            # proc = subprocess.call( cmd, shell=True, bufsize=0,) # COMMAND LINE
        print('success: bet_brains_fmap')
    
    def prepare_fmap(self, echo_tdiff=0.00246):
        """ Need to prepare field map before B0 unwarping
        https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide#SIEMENS_data
        Run after brain extraction on magnitude image
        Note too important which magnitude image is chosen
        Echo time difference need to be in ms (echotimes in JSON)
        fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE (in ms)> [--nocheck]
        The output of the Fsl_prepare_fieldmap GUI is a correctly calibrated fieldmap in units of rad/s     
        **
        echo_tdiff = 0.00246 # seconds, difference between echo times of field map magnitude images [0.0047,0.00716]
        **
        **Values to enter into FEAT**:
        EPI TE (s) of localizer = 0.0396
        'dwell_time' refers to effective echo spacing of EPI data (s) = 0.000580009
        NOTE: effective echo spacing and EPI TE (ms) refer to the fMRI EPI data! not the field map data - look in the JSON file for the bold runs
        """
    
        echo_tdiff = echo_tdiff * 1000 # needs to be in milliseconds
        phase_image = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_{}_run-03_fmap.nii.gz'.format(self.subject,self.session))
        mag_image = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_{}_run-01_fmap_brain.nii.gz'.format(self.subject,self.session))
        outFile = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_{}_acq-fmap.nii.gz'.format(self.subject,self.session))
        # bet inFile outFile
        cmd = 'fsl_prepare_fieldmap {} {} {} {} {} [--nocheck]'.format('SIEMENS',phase_image,mag_image,outFile,echo_tdiff)
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        print('success: prepare_fmap')
