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

import os, subprocess, sys, glob
import shutil as sh
import nibabel as nib
import pandas as pd
import numpy as np
from IPython import embed as shell # for Oly's debugging only

class preprocess_class(object):
    """
    Define a class for the preprocessing of the MRI data.
    
    Args:
        subject (str): Subject number main experiment.
        mri_subject (str): The chronological subject number from scanner.
        session (str): Current session out of two.
        analysis_dir (str): Path to the analysis scripts.
        source_dir (str): Path to the DICOM files for the MRI acquisition.
        raw_dir (str): Path to the NIFTI files after converstion from DICOM.
        deriv_dir (str): Path to the derivatives directory.
        mask_dir (str): Path to the directory with MRI masks.
        template_dir (str): Path to the directory with the FEAT/FSL template files.
        FWHM (int): Smoothing kernel for localizer analysis (full-width half maximum)
        t1_sess (str): Session to use for T1 in FEAT.
        
    Attributes:
        subject (str): Subject number main experiment.
        mri_subject (str): The chronological subject number from scanner.
        session (str): Current session out of two.
        analysis_dir (str): Path to the analysis scripts.
        source_dir (str): Path to the DICOM files for the MRI acquisition.
        raw_dir (str): Path to the NIFTI files after converstion from DICOM.
        deriv_dir (str): Path to the derivatives directory.
        mask_dir (str): Path to the directory with MRI masks.
        template_dir (str): Path to the directory with the FEAT/FSL template files.
        FWHM (str): Smoothing kernel for localizer analysis (full-width half maximum)
        t1_sess (str): Session to use for T1 in FEAT.
        preprocess_dir (str): Path to the preprocessing directory.
        preprocessing_job_path (str): Path to the text file for batch scripting (Unix).
    
    Methods:
        __init__(subject, mri_subject, session, analysis_dir, source_dir, raw_dir, deriv_dir, mask_dir, template_dir, FWHM, t1_path)
        dicom2bids()
        rename_logfiles(behav_sess)
        remove_timestamps_logfiles(behav_sess)
        copy_logfiles(behav_sess)
        housekeeping()
        raw_copy()
        bet_brains_T1( postfix='brain')
        preprocess_fsf()
        register_native2native_target(task='rsa', session='ses-01', bold_run='run-01')
    """
    
    def __init__(self, subject, mri_subject, session, analysis_dir, source_dir, raw_dir, deriv_dir, mask_dir, template_dir, EPI_TE, EPI_EECHO, TDIFF_ECHO, UNWARP, FWHM, t1_sess):        
        self.subject        = 'sub-'+str(subject)
        self.mri_subject    = 'sub-'+str(mri_subject)
        self.session        = str(session)
        self.analysis_dir   = str(analysis_dir)
        self.source_dir     = str(source_dir)   # dicoms
        self.raw_dir        = str(raw_dir)      # raw niftis
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.EPI_TE         = str(EPI_TE)
        self.EPI_EECHO      = str(EPI_EECHO)
        self.TDIFF_ECHO     = str(TDIFF_ECHO)
        self.UNWARP         = str(UNWARP)
        self.FWHM           = str(FWHM)
        self.t1_sess        = str(t1_sess) # with _brain extension
            
        if not os.path.isdir(self.raw_dir):
            os.mkdir(self.raw_dir)
        if not os.path.isdir(self.deriv_dir):
            os.mkdir(self.deriv_dir)
        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
        
        self.preprocess_dir = os.path.join(self.deriv_dir, 'preprocessing')
        if not os.path.isdir(self.preprocess_dir):
            os.mkdir(self.preprocess_dir)
        
        if not os.path.isdir(os.path.join(self.preprocess_dir, 'task-colors')):
            os.mkdir(os.path.join(self.preprocess_dir, 'task-colors'))
        if not os.path.isdir(os.path.join(self.preprocess_dir, 'task-letters')):
            os.mkdir(os.path.join(self.preprocess_dir, 'task-letters'))
        if not os.path.isdir(os.path.join(self.preprocess_dir, 'task-rsa')):
            os.mkdir(os.path.join(self.preprocess_dir, 'task-rsa'))
            
        # path to processed physiology files
        self.phys_dir = os.path.join(self.raw_dir, 'physiology_processed')
                
        # write unix commands to job to run in parallel
        self.preprocessing_job_path = os.path.join(self.analysis_dir, 'jobs', 'job_preprocessing_{}.txt'.format(self.subject))
        if self.session == 'ses-01':
            self.preprocessing_job = open(self.preprocessing_job_path, "w")
            self.preprocessing_job.write("#!/bin/bash\n")
            self.preprocessing_job.close()
    
    
    def dicom2bids(self, ):
        """Convert MRI data from dicom to nifti format and structure files according to the BIDS standard. 
        
        Notes:
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
        
        dicom_dir = os.path.join(self.source_dir, self.mri_subject, mri_session) #/Users/olympiacolizoli/Aeneas/mountpoint6/raw
        config_file  = os.path.join(self.analysis_dir, 'config', config)
        output_dir = os.path.join(self.raw_dir, 'nifti')

        # Create folder structure if not yet existing
        if not os.path.isdir(os.path.join(self.raw_dir, 'nifti')):
            os.mkdir(os.path.join(self.raw_dir, 'nifti'))

        cmd = 'dcm2bids -d {} -p {} -s {} -c {} -o {}'.format(dicom_dir, self.subject, self.session, config_file, output_dir)
        
        results = subprocess.call(cmd, shell=True, bufsize=0)
        print('converting {} convert2bids_letter-color_main_{}.json'.format(self.mri_subject,self.subject))
        print('success: dicom2bids')
    
    
    def rename_logfiles(self, behav_sess):
        """Rename the log files to match MRI files: e.g., sess-2 -> ses-02.
        
        Args:
            behav_sess (str): The session string within the behavioral logfile names.
        """
        ###################
        # RENAME 'behav' 'sess-1' to 'ses-01' in logfiles/behav folder before copying
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'behav')
        for f in os.listdir(dir_path):
            if behav_sess in f:
                f_new = f.replace(behav_sess, self.session)
                os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                print('old={} , new={}'.format(f,f_new))
        
        ###################
        # RENAME 'func' 'sess-1' to 'ses-01' in logfiles/func folder before copying
        # then runs 'run-1' to 'run-01' etc
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
        behav_runs = ['run-1','run-2','run-3','run-4']
        bids_runs = ['run-01','run-02','run-03','run-04']
        # loops through functional files
        for f in os.listdir(dir_path):
            # replace sessions
            if behav_sess in f:
                f_new = f.replace(behav_sess, self.session)
                os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                print('old={} , new={}'.format(f,f_new))
        # restart loop because new session changed names
        for f in os.listdir(dir_path):
            # replaces runs
            for r,this_run in enumerate(behav_runs):
                if this_run in f:
                    f_new = F.replace(this_run, bids_runs[r])
                    os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                    print('old={} , new={}'.format(f,f_new))
        
        ###################
        # RENAME localizer runs only: remove '_run-01' because only 1 run of each
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
        # loops through functional files
        for F in os.listdir(dir_path):
            # remove run from localizer event file names
            if ('task-colorsloc' in f) or ('task-lettersloc' in f):
                f_new = f.replace('_run-01', '')
                os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                print('old={} , new={}'.format(f,f_new))
                try: # some files have runs 5 or 6
                    f_new = f.replace('_run-5', '')
                    os.rename(os.path.join(dir_path,F),os.path.join(dir_path,f_new)) # old,new
                    f_new = f.replace('_run-6', '')
                    os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                except:
                    pass
        
        ################### 
        # RENAME localizer runs: 'task-colorsloc' to 'task-colors'
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
        # loops through functional files
        for f in os.listdir(dir_path):
            # remove run from localizer event file names
            if ('task-colorsloc' in f):
                f_new = f.replace('task-colorsloc', 'task-colors') # old, new
                os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                print('old={} , new={}'.format(f,f_new))
            elif ('task-lettersloc' in f):
                f_new = f.replace('task-lettersloc', 'task-letters') # old, new
                os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                print('old={} , new={}'.format(f,f_new))
                    
        print('success: rename_logfiles')
    
    
    def remove_timestamps_logfiles(self, behav_sess):
        """Remove the trailing timestamps from the psychopy output: last 16 characters _00200312-125726.
        
        Args:
            behav_sess (str): The session string within the behavioral logfile names.
        """
        T = 16 # length of trailing characters to remove
        
        ###################
        # behav folder
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'behav')
        # loops through functional files
        for f in os.listdir(dir_path):
            if ('task-iconic' in f) or ('prefs' in f) or ('task-consistency' in f):
                Fname, extension = os.path.splitext(f) # split name from extension
                F_new = Fname[:-T]+extension # remove T trailing T characters, add extension back
                os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
                print('old={} , new={}'.format(f,f_new))
                
        ###################
        # func folder
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
        # loops through functional files
        for f in os.listdir(dir_path):
            fname, extension = os.path.splitext(f) # split name from extension
            f_new = fname[:-T]+extension # remove T trailing T characters, add extension back
            os.rename(os.path.join(dir_path,f),os.path.join(dir_path,f_new)) # old,new
            print('old={} , new={}'.format(f,f_new))
        
        print('success: remove_timestamps_logfiles')
        
    def copy_logfiles(self, behav_sess):
        """Copy events and log files into the bids_raw directory with the NIFTI files.
        
        Args:
            behav_sess (str): The session string within the behavioral logfile names.
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
        for f in os.listdir(os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')):
            src = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func', f)
            dst = os.path.join(self.raw_dir, 'nifti', self.subject, self.session, 'func', f)
            sh.copy(src, dst)
        print('func: src={} , dst={}'.format(src,dst))
        
        print('success: copy_logfiles')
    
    def copy_physio(self, ):
        """Copy the physio files into the bids_raw directory.
        
        Notes: 
            First run the separate script process_source_physio.py
            Already renamed and runs are matched.
        """
        ###################
        # copy subject's files current session
        ###################
        
        subj_fns = glob.glob(os.path.join(self.phys_dir, '{}_{}_*_physio.tsv'.format(self.subject,self.session)))
        
        for phys in subj_fns:
            src = os.path.join(self.phys_dir, phys) # from physiology_processed
            dst = os.path.join(self.raw_dir, self.subject, self.session, 'func', phys.split("/")[-1]) # bids_raw 
            sh.copyfile(src, dst)
            print('copying: src={} , dst={}'.format(src,dst))
        print('success: copy_physio')
    
    
    def remove_filename(self, ):
        
        subj_dir = os.path.join(self.deriv_dir, self.subject, self.session, 'func')
        
        cmd = 'rm {}*physio.tsv'.format(subj_dir)
        print(cmd)
        #results = subprocess.call(cmd, shell=True, bufsize=0)
        print('success: remove_filename')
        
    
    def housekeeping(self, ):
        """Rename then copies events, log, and physio files into the bids_raw directory with the NIFTI files.
        
        Notes:
            WARNING!! don't remove timestamps twice on same participant by rerunning remove_timestamps_logfiles().
        """
        if self.session == 'ses-01':
            behav_sess = 'sess-1' #old 
        else:
            behav_sess = 'sess-2'
        
        # self.rename_logfiles(behav_sess)            # rename session and runs
        # ## WARNING!! don't remove timestamps twice on same participant
        # self.remove_timestamps_logfiles(behav_sess) # removes the trailing timestamps from the psychopy output
        # self.copy_logfiles(behav_sess)              # copy logfiles to bids_raw folder with mri data
        # self.copy_physio()                            # copy processed physio files into bids_raw directory after process_source_physio()
        print('success: housekeeping')
    
    
    def raw_copy(self,):
        """Copy the entire bids_raw directory into derivatives folder for further processing/analysis.
        
        Notes:
            The bids_raw directory should be backed up on the DAC. 
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


    def sbref_raw_copy(self,):
        """Copy the sbref images from the bids_raw directory into derivatives folder for further processing/analysis.
        
        Notes:
            The bids_raw directory should be backed up on the DAC. 
            All further analysis should be run within the derivatives folder.
        """
        
        subj_fns = glob.glob(os.path.join(self.raw_dir, 'nifti', self.subject, self.session, 'func', '{}_{}_*'.format(self.subject,self.session)))
        
        for sbref in subj_fns:
            src = os.path.join(self.phys_dir, sbref) # from physiology_processed
            dst = os.path.join(self.deriv_dir, self.subject, self.session, 'func', sbref.split("/")[-1]) # bids_raw 
            sh.copyfile(src, dst)
        print('success: sbref_raw_copy')
        
        
    def bet_brains_T1(self, postfix='brain'):
        """Run brain extraction on all T1s.
        
        Args:
            postfix (str): The string to add to the end of the brain-extracted images (default 'brain').
        
        Notes:
            Always check visually in fsleyes.
        """
        mri_in = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w.nii.gz'.format(self.subject, self.session))
        
        if os.path.exists(inFile):
            mri_out = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w_{}.nii.gz'.format(self.subject, self.session, postfix))
            # bet inFile outFile
            cmd = 'bet {} {}'.format(mri_in, mri_out)
            print(cmd)
            results = subprocess.call(cmd, shell=True, bufsize=0)

            # reorient to mni space if necessary (sometimes turned around)
            # cmd = 'fslreorient2std {} {}'.format(outFile,outFile)
            # proc = subprocess.call( cmd, shell=True, bufsize=0,) # COMMAND LINE
        else:
            print('No T1 for {} {}'.format(self.subject, self.session))
        print('success: bet_brains_T1')
    
    
    def bet_brains_fmap(self, postfix='brain'):
        """Run brain extraction on magnitude fieldmap (fmap) images.
        
        Args:
            postfix (str): The string to add to the end of the brain-extracted images (default 'brain').
        
        Notes:
            Important! Always check visually in fsleyes.
            Brain extraction needs to be 'tight', meaning it is better to exclude brain voxels than include noisy non-brain voxels.
        """
        if int(self.UNWARP):
            for run in ['01','02','03']: # first 2 images are magnitude, 3rd image is phase
                inFile = os.path.join(self.deriv_dir,  self.subject,  self.session, 'fmap', '{}_{}_run-{}_fmap.nii.gz'.format(self.subject,self.session,run))
                outFile = os.path.join(self.deriv_dir,  self.subject,  self.session, 'fmap', '{}_{}_run-{}_fmap_{}.nii.gz'.format(self.subject,self.session,run,postfix))
                # bet inFile outFile
                cmd = 'bet {} {} -f 0.6'.format(inFile, outFile) 
                print(cmd)
                results = subprocess.call(cmd, shell=True, bufsize=0)
    
                # reorient to mni space if necessary (sometimes turned around)
                # cmd = 'fslreorient2std {} {}'.format(outFile,outFile)
                # proc = subprocess.call( cmd, shell=True, bufsize=0,) # COMMAND LINE
        else:
            print('No fieldmap for {} {}'.format(self.subject, self.session))
        print('success: bet_brains_fmap')
    
    
    def prepare_fmap(self, ):
        """Prepare field map before B0 unwarping.
        
        Notes:
            https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide#SIEMENS_data
            Run after brain extraction on magnitude image.
            Not too important which magnitude image is chosen.
            Echo time difference need to be in ms (echotimes in JSON).
            fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE (in ms)> [--nocheck]
            The output of the Fsl_prepare_fieldmap GUI is a correctly calibrated fieldmap in units of rad/s     
            **
            echo_tdiff = 0.00246 # seconds, difference between echo times of field map magnitude images [0.0047,0.00716]
            **
            **
            Values to enter into FEAT during unwarping step in preprocessing:
            **
            EPI TE (s) of localizer = 0.0396
            'dwell_time' refers to effective echo spacing of EPI data (s) = 0.000580009
            **
            Important! effective echo spacing and EPI TE (ms) refer to the fMRI EPI data! not the field map data - look in the JSON file for the bold runs.
            No in-plane acceleration was used (this is different from multi-band acceleration.). In plane acceleration (GRAPPA) would be indicated in the 
            PDF of the Siemens sequence from the Donders Intranet under "Resolution - iPAT". 
            https://intranet.donders.ru.nl/fileadmin/user_upload/DCCN/Laboratories/MRI/MR_protocols/Prisma/cmrr_2p0iso_mb4_TR1500.pdf
        """
        if int(self.UNWARP):
            phase_image = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_{}_run-03_fmap.nii.gz'.format(self.subject,self.session))
            mag_image = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_{}_run-01_fmap_brain.nii.gz'.format(self.subject,self.session))
            outFile = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_{}_acq-fmap.nii.gz'.format(self.subject,self.session))
            # bet inFile outFile
            cmd = 'fsl_prepare_fieldmap {} {} {} {} {} [--nocheck]'.format('SIEMENS', phase_image, mag_image, outFile, self.TDIFF_ECHO)
            print(cmd)
            results = subprocess.call(cmd, shell=True, bufsize=0)
        else:
            print('No fieldmap for {} {}'.format(self.subject, self.session))
        print('success: prepare_fmap')


    def preprocess_fsf(self, native_target_session='ses-01', native_target_task_run='rsa_run-01'):
        """Create the FSF files for each subject's preprocessing in FEAT: 
        B0 unwarping, nonlinear reg to MNI, motion correction, high pass filtering, brain extraction, smoothing (localizers).
        
        Args:
            native_target_session (str): the session used for creating the native target image.
            native_target_task_run (str): the task and run used for creating the native target image.
        
        Notes:
            Smoothing is on only for the localizers.
            Registration to standard space is only on for the native_target.
            Uses the preprocessing_template.fsf in analysis/templates to search and replace the markers.
        """    
        
        template_filename = os.path.join(self.template_dir,'preprocessing_template.fsf')
        
        markers = [
            '[$OUTPUT_PATH]', 
            '[$NR_TRS]',        # number of volumes
            '[$NR_VOXELS]',     # total number of voxels
            '[$FWHM]',
            '[$UNWARP_ON]',     # turn on unwarping
            '[$EPI_EECHO]',     # EPI DWELL TIME, EFFECTIVE ECHO SPACING
            '[$EPI_TE]',        # EPI echo time
            '[$INPUT_FILENAME]', # BOLD data
            '[$FMAP]',
            '[$FMAP_MAG_BRAIN]',
            '[$SBREF_BRAIN]',
            '[$T1_BRAIN]',
            '[$MNI_BRAIN]',
            '[$REG_MNI_BRAIN]'  # on only for native target
        ]
        
        # all nii.gz files in func folder need to be preprocessed in FEAT
        for task_run in ['letters','colors','rsa_run-01','rsa_run-02','rsa_run-03','rsa_run-04']:
            
            if 'rsa' in task_run:
                task = task_run.split("_")[0]
            else:
                task = task_run
                
            bold_in = os.path.join(self.deriv_dir, self.subject, self.session, 'func', '{}_{}_task-{}_bold.nii.gz'.format(self.subject,self.session,task_run))

            if os.path.exists(bold_in):  ##### SKIP IF BOLD FILE DOES NOT EXIST #####

                FSF_filename    = os.path.join(self.preprocess_dir,'task-{}'.format(task), '{}_{}_task-{}_preprocessing.fsf'.format(self.subject, self.session, task_run) ) # save fsf
                output_path     = os.path.join(self.preprocess_dir,'task-{}'.format(task), '{}_{}_task-{}_preprocessing'.format(self.subject, self.session, task_run) ) 
                
                nii = nib.load(bold_in).get_data() # only do once 
                nr_trs = str(nii.shape[-1])
                nr_voxels = str(nii.size)

                # replacements
                sbref_brain     = os.path.join(self.deriv_dir,self.subject,self.session,'func','{}_{}_task-{}_bold_sbref.nii.gz'.format(self.subject, self.session, task_run))
                fmap            = os.path.join(self.deriv_dir,self.subject,self.session,'fmap','{}_{}_acq-fmap.nii.gz'.format(self.subject,self.session))
                fmap_mag_brain  = os.path.join(self.deriv_dir,self.subject,self.session,'fmap','{}_{}_run-01_fmap_brain.nii.gz'.format(self.subject,self.session))
                t1_brain        = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w_brain.nii.gz'.format(self.subject, self.t1_sess))
                mni_brain       = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz')
                
                # conditional preprocessing
                if task == 'rsa':
                    FWHM = '0' # turn smoothing off for mulitvariate analyses
                else:
                    FWHM = self.FWHM
                
                # registration to standard space only for native target
                if (native_target_task_run in task_run) and (native_target_session in self.session):
                    REG_MNI_BRAIN = '1' 
                else:
                    REG_MNI_BRAIN = '0'
                    
                replacements = [ # needs to match order of 'markers'
                    output_path,
                    nr_trs,
                    nr_voxels,
                    FWHM,
                    self.UNWARP,# true=ON
                    self.EPI_EECHO, # dwell time is effective echo spacing (EPI data not field map!!)
                    self.EPI_TE,
                    bold_in,
                    fmap,
                    fmap_mag_brain,
                    sbref_brain,
                    t1_brain,
                    mni_brain,
                    REG_MNI_BRAIN]

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
                self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                self.preprocessing_job.write(cmd)   # feat command
                self.preprocessing_job.write("\n\n")  # new line
                self.preprocessing_job.close()
                print('success: {}'.format(FSF_filename))
            else:
                print('cannot make FSF: missing {}'.format(BOLD))
        
    
    def register_native2native_target(self, task='rsa', session='ses-01', bold_run='run-01'):
        """Calcuate the linear transformations for individual runs into native target space. 
        
        Args:
            task (str): The task to use for the registration target ('rsa' or 'loc').
            session (str): The session to use for the registration target ('ses-01' or 'ses-02').
            bold_run (str): The bold run to use for registration target.
        
        Notes:
        ------
            Only compute, not apply the transforms!
            The target should be the first session, 1st run (RSA1), target of motion correction, middle volume ("example_func").
        """        
        
        for preprocess_task in ['letters','colors','rsa']:
            
            dir_path = os.path.join(self.preprocess_dir, 'task-{}'.format(preprocess_task) ) # preprocessing directory path
            
            for f in os.listdir(dir_path):
                if (self.session in f) and ('.feat' in f): # check if feat directory and avoid double sessions
                    
                    # skip if already native_target run
                    if not ((task in preprocess_task) and (session in f) and (bold_run in f)):
                        
                        mri_in        = os.path.join(dir_path, f, 'example_func.nii.gz')
                        mri_out       = os.path.join(dir_path, f, 'reg', 'example_func2native_target')
                        mri_inverse   = os.path.join(dir_path, f, 'reg', 'native_target2example_func')
                        native_target = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'example_func.nii.gz')
        
                
                        ###################
                        # FLIRT COMMAND
                        ###################
                        # flirt -in {} -ref {} -out {}.nii.gz -omat {}.mat -bins 256 -cost mutualinfo -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp sinc -sincwidth 7 -sincwindow hanning
                        cmd1 = 'flirt -in {} -ref {} -out {}.nii.gz -omat {}.mat -bins 256 -cost mutualinfo -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp sinc -sincwidth 7 -sincwindow hanning'.format(mri_in, native_target, mri_out, mri_out)
                        print(cmd1)
                        # results = subprocess.call(cmd, shell=True, bufsize=0)
                        
                        ###################
                        # INVERT TRANSFORM
                        ###################
                        # convert_xfm -omat {}.mat -inverse {}.mat
                        cmd2 = 'convert_xfm -omat {}.mat -inverse {}.mat'.format(mri_inverse, mri_out)
                        print(cmd2)
                        
                        # open preprocessing job and write commands as new line
                        self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                        self.preprocessing_job.write(cmd1)   # command
                        self.preprocessing_job.write("\n\n")  # new line
                        self.preprocessing_job.write(cmd2)   # command
                        self.preprocessing_job.write("\n\n")  # new line
                        self.preprocessing_job.close()                
        print('success: register_native2native_target')
    
    
    def invert_registrations(self, task='rsa', session='ses-01', bold_run='run-01'):
        """Invert the native_target to MNI non-linear registration and the native to native target registration.
        
        Args:
            task (str): The task to use for the registration target ('rsa' or 'loc').
            session (str): The session to use for the registration target ('ses-01' or 'ses-02').
            bold_run (str): The bold run to use for registration target.
        
        Notes:
        ------
            Only compute, not apply the transforms!
            The target should be the first session, 1st run (RSA1), target of motion correction, middle volume ("example_func").
            https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/UserGuide#A--inwarp
        """
        
        native_target               = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'example_func.nii.gz')
        native_target2standard_warp = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'reg', 'example_func2standard_warp.nii.gz')
        standard2native_target_warp = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'reg', 'standard2example_func_warp.nii.gz')
        
        ###################
        # INVERT FNIRT TRANSFORM NATIVE TARGET TO MNI
        ###################
        
        # invwarp --ref=my_struct --warp=warps_into_MNI_space --out=warps_into_my_struct_space
        cmd1 = 'invwarp --ref={} --warp={} --out={}'.format(native_target, native_target2standard_warp, standard2native_target_warp)
        print(cmd1)
                
        # open preprocessing job and write commands as new line
        self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
        self.preprocessing_job.write(cmd1)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.close()
        
        ###################
        # INVERT TRANSFORM example_func2native_target
        ###################
        for preprocess_task in ['letters','colors','rsa']:
            
            dir_path = os.path.join(self.preprocess_dir, 'task-{}'.format(preprocess_task) ) # preprocessing directory path
            
            for f in os.listdir(dir_path):
                if (self.session in f) and ('.feat' in f): # check if feat directory and avoid double sessions
                    
                    # skip if already native_target run
                    if not ((task in preprocess_task) and (session in f) and (bold_run in f)):
                        
                        example_func2native_target  = os.path.join(dir_path, f, 'reg', 'example_func2native_target')
                        native_target2native        = os.path.join(dir_path, f, 'reg', 'native_target2example_func') # native_target2native
                        
                        # CONVERT XFM
                        # convert_xfm -omat {}.mat -inverse {}.mat
                        cmd2 = 'convert_xfm -omat {}.mat -inverse {}.mat'.format(native_target2native, example_func2native_target)
                        print(cmd2)

                        # open preprocessing job and write commands as new line
                        self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                        self.preprocessing_job.write(cmd2)   # command
                        self.preprocessing_job.write("\n\n")  # new line
                        self.preprocessing_job.close()
        print('success: invert_registrations')
        
    
    def register_ventricle2native(self, task='rsa', session='ses-01', bold_run='run-01'):
        """Register the 4th ventricle (MNI space) to native space.
        
        Args:
            task (str): The task to use for the registration target ('rsa' or 'loc').
            session (str): The session to use for the registration target ('ses-01' or 'ses-02').
            bold_run (str): The bold run to use for registration target.
        
        Notes:
        ------
            Apply the transforms (first inverse warp, the inverse FLIRT transforms)
            The target should be the first session, 1st run (RSA1), target of motion correction, middle volume ("example_func").
            Threshold and mask ventricle (0.5 to maintain save number of voxels in native target)
        """
        
        native_target               = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'example_func.nii.gz')
        standard2native_target_warp = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'reg', 'standard2example_func_warp.nii.gz')
        ventricle                   = os.path.join(self.mask_dir, '4th_ventricle_MNI.nii.gz')
        ventricle2native_target     = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'reg', 'ventricle2example_func.nii.gz')
        
        ###################
        # APPLY INVERSE WARP 4th ventricle from MNI to native target (apply inverse warp)
        ###################
        # applywarp --ref=my_struct --in=ACC_left --warp=warps_into_my_struct_space --out=ACC_left_in_my_struct_space --interp=nn
        cmd1 = 'applywarp --ref={} --in={} --warp={} --out={} --interp=nn'.format(native_target, ventricle, standard2native_target_warp, ventricle2native_target)
        print(cmd1)
                
        # open preprocessing job and write commands as new line
        self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
        self.preprocessing_job.write(cmd1)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.close()
        
        ###################
        # APPLY INVERSE FLIRT 
        ###################
        for preprocess_task in ['letters','colors','rsa']:
            
            dir_path = os.path.join(self.preprocess_dir, 'task-{}'.format(preprocess_task) ) # preprocessing directory path
            
            for f in os.listdir(dir_path):
                if (self.session in f) and ('.feat' in f): # check if feat directory and avoid double sessions
                    
                    # skip if already native_target run
                    if not ((task in preprocess_task) and (session in f) and (bold_run in f)):
                        
                        native_target2native    = os.path.join(dir_path, f, 'reg', 'native_target2example_func') # native_target2native
                        ventricle2example_func  = os.path.join(dir_path, f, 'reg', 'ventricle2example_func.nii.gz')
                        example_func            = os.path.join(dir_path, f, 'example_func.nii.gz')
                        
                        # APPLY FLIRT
                        #flirt -in {} -applyxfm -init {} -out {} -paddingsize 0.0 -interp sinc -sincwidth 7 -sincwindow hanning -ref {}
                        
                        cmd2 = 'flirt -in {} -applyxfm -init {}.mat -out {} -paddingsize 0.0 -interp sinc -sincwidth 7 -sincwindow hanning -ref {}'.format(ventricle2native_target, native_target2native, ventricle2example_func, example_func)
                        print(cmd2)
                        
                        # THRESHOLD AND MASK VENTRICLE 
                        # fslmaths ventricle2example_func.nii.gz -thr 0.5 -bin ventricle2example_func.nii.gz 
                        cmd3 = 'fslmaths {} -thr 0.5 -bin {} '.format(ventricle2example_func, ventricle2example_func)
                        print(cmd3)
                        
                        # open preprocessing job and write commands as new line
                        self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                        self.preprocessing_job.write(cmd2)   # command
                        self.preprocessing_job.write("\n\n")  # new line
                        self.preprocessing_job.write(cmd3)   # command
                        self.preprocessing_job.write("\n\n")  # new line
                        self.preprocessing_job.close()
        print('success: register_ventricle2native')
       
        
    def ventricle_fluctations(self,dtype,session):
        # extracts the 4th ventricle fluctuations for the nuisance regression    
        # files have been copied to the combined folder
        BOLD_data = np.array(nib.load(os.path.join(self.project.base_dir, self.subject.initials, 'FSL_files', 'detection_data_{}_{}_{}.nii.gz'.format(dtype,self.subject.initials,session))).get_data(), dtype=float)
        mask = np.array(nib.load(os.path.join(self.project.base_dir, self.subject.initials, 'FSL_files', 'masks', '{}_{}.nii.gz'.format('4th_ventricle_MNI',session))).get_data(), dtype=bool)
        roi_data = bn.nanmean(BOLD_data[mask,:], axis=0)
        # roi_data = (roi_data - np.mean(roi_data)) / np.std(roi_data) # center variable around 0 mean
        return roi_data


    def physiological_noise_removal(self, ):
        """Remove physiological noise and 4th ventricle signal.
                
        Notes:
        ------
            https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PNM
            'Pnm_gui' will open the PNM GUI in FSL.
            commands are 'popp' then 'pnm_evs'
            Columns: 'idx', 'HR', 'Resp', 'trigger' but NO column headers!
            Sampling rate = 5000 Hz, TR = 1.5, Slice order = custom b/c mutli-band 4 interleaved up (Siemens), Scanner Slice Direction = z, 
            Orders: Cardiac=4, Resp=4, Cardiac_Int=2, Resp_Int=2
            Options: RVT, HeartRate, Smoothing defaults 
        """
        
        slice_timing = os.path.join(self.template_dir, 'bold_slice_timing.txt') # custom multi-band 4 interleaved up
        
        for rsa_run in ['run-01', 'run-02', 'run-03', 'run-04']:
            
            physio_input = os.path.join(self.deriv_dir, self.subject, self.session, 'func', '{}_{}_task-rsa_{}_physio.txt'.format(self.subject, self.session, rsa_run))
            
            if os.path.exists(phyio_input):
                # convert physio file to standard file format and name
                bold_input = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'filtered_func_data.nii.gz')
                pnm_input = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'pnm_input.txt')
                pnm_dir = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run))
                pnm_out = os.path.join(self.preprocess_dir, 'task-{}'.format(task), '{}_{}_task-{}_{}_preprocessing.feat'.format(self.subject, session, task, bold_run), 'pnm')
                
                ###### CONVERT TEXT FILE FORMAT
                # fslFixText sub-201_ses-01_task-rsa_run-01_physio.txt pnm_input.txt
                cmd1 = 'fslFixText {} {}'.format(this_physio, pnm_input)
                print(cmd1)
        
                ###### PROCESS PHYSIO
                # popp -i pnm_input.txt -o pnm -s 5000 --tr=1.5 --smoothcard=0.1 --smoothresp=0.1 --resp=4 --cardiac=3 --trigger=5 --rvt --heartrate --verbose
                cmd2 = 'popp -i {} -o pnm -s 5000 --tr=1.5 --smoothcard=0.1 --smoothresp=0.1 --resp=4 --cardiac=3 --trigger=5 --rvt --heartrate --verbose'.format(pnm_input)
                print(cmd2)
                
                ###### MAKE RETROICOR REGRESSORS (NIFTI)
                # pnm_evs --oc=4 --or=4 --multc=2 --multr=2 --slicedir='z' --slicetiming=bold_slice_timing.txt -v --tr=1.5 -i filtered_func_data.nii.gz -o pnmevs -r pnm_resp.txt -c pnm_card.txt --rvt=pnm_rvt.txt --heartrate=pnm_hr.txt --verbose
                cmd3 = 'pnm_evs --oc=4 --or=4 --multc=2 --multr=2 --slicedir=z --slicetiming={} --tr=1.5 -i {} -o {}/pnmevs -r {}/pnm_resp.txt -c {}/pnm_card.txt --rvt={}/pnm_rvt.txt --heartrate={}/pnm_hr.txt --verbose'.format(slice_timing, bold_input, pnm_base, pnm_base, pnm_base, pnm_base, pnm_base)
                print(cmd3)
                
                # open preprocessing job and write commands as new line
                self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                self.preprocessing_job.write(cmd1)   # command
                self.preprocessing_job.write("\n\n")  # new line
                self.preprocessing_job.write(cmd2)   # command
                self.preprocessing_job.write("\n\n")  # new line
                self.preprocessing_job.write(cmd3)   # command
                self.preprocessing_job.write("\n\n")  # new line
                self.preprocessing_job.close()
        
                # # grab regressors:
                # regressors = [reg for reg in np.sort(glob.glob(base + 'ev*.nii*'))]
                # text_file = open(base+'_evs_list.txt', 'w')
                # for reg in regressors:
                #     text_file.write('{}\n'.format(reg))
                # text_file.close()
        print('success: physiological_noise_removal')
        
        
    def GLM_nuisance(self,tsnr):
        # Nuisance GLM for the tSNR analysis:
        # for tSNR:
        # removing  i) the physiological fluctuations, ii) the task-evoked fluctuations, and iii) the 4th ventricle fluctuations
        # for task-evoked 
        # removing  i) the physiological fluctuations, ii) the 4th ventricle fluctuations
        # Include all of the above regressors in one nuisance GLM (i.e. in one step along with the RETROICOR). 
        # run from combined folder (except brain masks)

        session = [self.runList[i] for i in self.conditionDict['task']][0].session
        dtype = 'hp'
        if tsnr:
            task_regressors = self.make_task_evoked_regressors(session) # make 3 column format from stimulus and feedback onsets
            dtype_out = 'hp'
           # dtype = 'hp'
        else:
            dtype_out = 'hp_evoked'
           # dtype = 'psc'
        ventricle_regressor = self.ventricle_fluctations(dtype=dtype,session=session)  # one column of average activity in ventricle per TR

        nr_runs = len([self.runList[i] for i in self.conditionDict['task']])
        nii_file = nib.load(os.path.join(self.project.base_dir, self.subject.initials, 'FSL_files', 'detection_data_{}_{}_{}.nii.gz'.format(dtype,self.subject.initials, session)))
        nii_file_data = nii_file.get_data()

        if tsnr:
            # make box car for task regressors, convolve with HRF
            task_regressors = self.convolve_hrf(task_regressors, nii_file_data.shape[-1], self.tr)
            # fig = plt.figure(figsize=(5,5))
            # ax = fig.add_subplot(211)
            # ax.plot(task_regressors[0][:400])
            # ax = fig.add_subplot(212)
            # ax.plot(task_regressors[1][:400])
            # fig.savefig(os.path.join(self.project.base_dir, self.subject.initials, 'FSL_files','task_regressors_{}.pdf'.format(session)))
        with open(os.path.join(self.project.base_dir, self.subject.initials, 'FSL_files', 'retroicor_data_{}_{}_evs_list.txt'.format(self.subject.initials, session)), 'r') as f:
            retroicor_files = [x.strip('\n') for x in f.readlines()]
        retroicor_evs = []
        for r in retroicor_files:
            retroicor_evs.append(nib.load(r))

        brain_mask = np.array(nib.load(self.runFile(stage = 'processed/mri/masks/anat', base='brain_mask')).get_data(), dtype=bool)
        brain_mask[:,:,:] = True

        residuals = np.zeros(nii_file_data.shape)
        nr_slices = nii_file.shape[2]

        for s in range(nr_slices):
   
            slice_indices = (np.arange(nr_slices) == s)
            mask = brain_mask[:,:,slice_indices][:,:,0] # voxels x,y inplane
   
            retroicor_slice_wise_evs = []
           for r in retroicor_evs:
               retroicor_slice_wise_evs.append(r.get_data()[0,0,s,:])
           retroicor_slice_wise_evs = np.vstack(retroicor_slice_wise_evs) # design matrix for retriocor
       
           # combine all regressors
           joined_design_matrix = np.mat(retroicor_slice_wise_evs).T  # timepoints x regressors
           if tsnr: # include task-evoked fluctuations as nuisance
               joined_design_matrix = np.mat(np.vstack([task_regressors, retroicor_slice_wise_evs, ventricle_regressor]).T)
           else:
               joined_design_matrix = np.mat(np.vstack([retroicor_slice_wise_evs, ventricle_regressor]).T) # X = time x regressors
           data = nii_file_data[mask,s,:].T  # time x voxels (flattened inplane x*y)
   
           # T is transpose, I inverse, # betas = # regressors
           betas = ((joined_design_matrix.T * joined_design_matrix).I * joined_design_matrix.T) * np.mat(data.T).T # betas = regressors x voxels
           explained_signal = np.array((np.mat(joined_design_matrix) * np.mat(betas))) # explained = time x voxels
           res = data - explained_signal # res = time x voxels

           # add mean over time back into residuals
           for tr in range(res.shape[0]): # time
               res[tr,:] = res[tr,:] + np.mean(data,axis=0)
           residuals[mask,s,:] = res.T # voxels x time
           # residuals shape is x-voxels x y-voxels x slices x time
   
        res_nii_file = nib.Nifti1Image(residuals, affine=nii_file.get_affine(), header=nii_file.get_header())
        res_nii_file.set_data_dtype(np.float32)
        nib.save(res_nii_file, os.path.join(self.project.base_dir, self.subject.initials, 'FSL_files', 'detection_data_{}_clean_{}_{}.nii.gz'.format(dtype_out,self.subject.initials, session)))     




    #### NOT USING

    
    def native_target_2_mni(self, ):
        """Compute the registration from the subject-specific native target to MNI space via the T1.
        
        Notes:
            FLIRT to T1 then FNIRT to MNI. Apply warp to native_target.
            Save transformation matrices to applywarp later.
            See: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/UserGuide#Example_uses
            >> flirt -ref my_betted_structural -in my_functional -dof 6 -omat func2struct.mat
            >> flirt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -in my_betted_structural -omat my_affine_transf.mat
            >> fnirt --in=my_structural --aff=my_affine_transf.mat --cout=my_nonlinear_transf --config=T1_2_MNI152_2mm
            >> applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_2mm --in=my_functional --warp=my_nonlinear_transf --premat=func2struct.mat --out=my_warped_functional
        """
        example_func = self.native_target + '_brain'
        sbref = os.path.join(self.deriv_dir,self.subject,self.session,'func','{}_{}_task-{}_{}_bold_sbref.nii.gz'.format(self.subject, self.session, 'rsa', 'run-03'))
        t1         = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w.nii.gz'.format(self.subject, self.t1_sess))
        t1_brain   = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w_brain.nii.gz'.format(self.subject, self.t1_sess))
        mni        = os.path.join(self.mask_dir, 'MNI152_T1_2mm.nii.gz') 
        mni_brain  = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz') 
    
        # define paths and output filenames
        example_func2sbref    = os.path.join(self.native_target_dir, 'example_func2sbref') 
        sbref2highres         = os.path.join(self.native_target_dir, 'sbref2highres') 
        
        example_func2highres    = os.path.join(self.native_target_dir, 'example_func2highres') 
        highres2mni_affine      = os.path.join(self.native_target_dir, 'highres2mni_affine') 
        highres2mni_fnirt       = os.path.join(self.native_target_dir, 'highres2mni_fnirt') 
        warpvol                 = os.path.join(self.native_target_dir, 'warpvol')
        example_func2mni        = os.path.join(self.native_target_dir, 'example_func2mni')
        
        # EPI target to SBREF
        cmd1 = 'flirt -ref {} -in {}.nii.gz -omat {}.mat -out {}.nii.gz -cost mutualinfo -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear'.format(sbref, example_func, example_func2sbref, example_func2sbref)
        print(cmd1)
        
        # SBREF to T1
        cmd2 = 'flirt -ref {} -in {} -omat {}.mat -out {}.nii.gz -cost mutualinfo -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear'.format(t1_brain, sbref, sbref2highres, sbref2highres)
        print(cmd2)
        
        # Concatenatve EPI target -> SBREF -> T1: to make EPI target to T1
        cmd3 = 'convert_xfm -omat {}.mat -concat {}.mat {}.mat'.format(example_func2highres, sbref2highres, example_func2sbref)
        print(cmd3)
        
        # Apply transform EPI target to T1
        cmd4 = 'flirt -in {}.nii.gz -applyxfm -init {}.mat -out {}.nii.gz -paddingsize 0.0 -interp trilinear -ref {}'.format(example_func, example_func2highres, example_func2highres, t1)
        print(cmd4)
        
        
        # EPI target to T1
        # cmd4 = 'flirt -ref {} -in {}.nii.gz -omat {}.mat -out {}.nii.gz -cost mutualinfo -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear'.format(t1_brain, example_func, example_func2highres, example_func2highres)
        # print(cmd1)
        # results = subprocess.call(cmd1, shell=True, bufsize=0)
    
        # T1 -> MNI
        cmd5 = 'flirt -ref {} -in {} -omat {}.mat -out {}.nii.gz'.format(mni_brain, t1_brain, highres2mni_affine, highres2mni_affine)
        print(cmd5)
        # results = subprocess.call(cmd2, shell=True, bufsize=0)
    
        # NON LINEAR
        # T1 > MNI space -> create FNIRT warpfile
        cmd6 = 'fnirt --ref{} --in={} --aff={}.mat --cout={}.nii.gz --iout={}.nii.gz'.format(mni, t1, highres2mni_affine, warpvol, highres2mni_fnirt)
        print(cmd6)
        # results = subprocess.call(cmd3, shell=True, bufsize=0)
    
        # NON LINEAR
        # EPI time series (MNI space) -> apply FNIRT warpfile
        cmd7 = 'applywarp --ref={} --in={}.nii.gz --warp={}.nii.gz --premat={}.mat --out={}.nii.gz'.format(mni, example_func, warpvol, example_func2highres, example_func2mni)
        print(cmd7)
        # results = subprocess.call(cmd4, shell=True, bufsize=0)
                    
        # open preprocessing job and write commands as new line
        self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
        self.preprocessing_job.write(cmd1)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.write(cmd2)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.write(cmd3)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.write(cmd4)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.write(cmd5)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.write(cmd6)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.write(cmd7)   # command
        self.preprocessing_job.write("\n\n")  # new line
        self.preprocessing_job.close()
        print('success: native_target_2_mni')
        
        
    def transform_2_mni(self, ):
        """Apply the same registration of native target-> MNI space to transform the filtered_func_data.nii.gz to MNI space (non-linear)
        
        Notes:
            FLIRT and FNIRT registrations have already been calculated based on native_target.
            Here, we are just applying the existing FNIRT warps to the preprocessed data
            See here: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/FAQ#How_do_I_transform_a_mask_with_FLIRT_from_one_space_to_another.3F
            https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/UserGuide#Now_what.3F_--_applywarp.21
        """
        
        # all nii.gz files in func folder need to be preprocessed in FEAT
        for task in ['letters','colors','rsa_run-01','rsa_run-02','rsa_run-03','rsa_run-04']:
            ##### SKIP IF BOLD INPUT FILE DOES NOT EXIST #####
            check_bold = os.path.join(self.deriv_dir, self.subject, self.session, 'func', '{}_{}_task-{}{}_bold.nii.gz'.format(self.subject,self.session,task,bold_run))
            
            if os.path.exists(check_bold):  ##### SKIP IF BOLD FILE DOES NOT EXIST #####

                epi     = os.path.join(self.preprocess_dir, 'task-{}'.format(task),'task-{}_{}_{}{}.feat'.format(task,self.subject,self.session,bold_run),'filtered_func_data.nii.gz')
                epi_mni = os.path.join(self.preprocess_dir, 'task-{}'.format(task),'task-{}_{}_{}{}.feat'.format(task,self.subject,self.session,bold_run),'filtered_func_data_mni.nii.gz')
                mni     = os.path.join(self.mask_dir, 'MNI152_T1_2mm.nii.gz') 
                example_func2highres    = os.path.join(self.native_target_dir, 'example_func2highres') 
                warpvol                 = os.path.join(self.native_target_dir, 'warpfile')

                # Apply FNIRT warpfile (non-linear)
                # EPI time series (in MNI space) -> apply FNIRT warpfile based on T1 (in MNI space) -> warped to MNI space
                cmd = 'applywarp --ref={} --in={}.nii.gz --warp={}.nii.gz --premat={}.mat --out={}.nii.gz'.format(mni, epi, warpvol, example_func2highres, epi_mni)
                print(cmd)
                # results = subprocess.call(cmd, shell=True, bufsize=0)
                # open preprocessing job and write command as new line
                self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                self.preprocessing_job.write(cmd)   # command
                self.preprocessing_job.write("\n\n")  # new line
                self.preprocessing_job.close()
        print('success: transform_2_mni')




    def motion_correction(self, ):
        """Run motion correction with subject-specific native target as reference volume.
        
        Notes:
            Call MCFLIRT from command line (FEAT has no option for overriding default reference volume.)
            Save motion matrices to add to first level GLM.
        """
        # go into subject's session's func folder and motion correct all nii.gz files
        for bold in os.listdir(os.path.join(self.deriv_dir,self.subject,self.session,'func')):
            if '.nii.gz' in bold:          
                # save in preprocessing folder of task, have to split 'nii.gz' twice to get base
                fname, extension = os.path.splitext(bold)
                fname, extension = os.path.splitext(fname)
                bold_out = os.path.join(self.motion_dir, fname+'_mcf')
                bold_in = os.path.join(self.deriv_dir,self.subject,self.session,'func', bold)
                ###################
                # motion correction (reference = native target)
                ###################
                cmd = 'mcflirt -in {} -r {}.nii.gz -stages 4 -sinc_final -o {} -stats -mats -plots'.format(bold_in, self.native_target, bold_out)
                print(cmd)
                # results = subprocess.call(cmd, shell=True, bufsize=0)

                # open preprocessing job and write command as new line
                self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                self.preprocessing_job.write(cmd)   # command
                self.preprocessing_job.write("\n\n")  # new line
                self.preprocessing_job.close()                
        print('success: motion_correction')