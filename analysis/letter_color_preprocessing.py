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
    def __init__(self, subject, mri_subject, session, analysis_dir, source_dir, raw_dir, deriv_dir, mask_dir, template_dir, timing_files_dir, EPI_TE, EPI_EECHO, FWHM):        
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
        self.EPI_TE         = str(EPI_TE)
        self.EPI_EECHO      = str(EPI_EECHO)
        self.FWHM           = str(FWHM)
            
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
        
        self.preprocess_dir = os.path.join(self.deriv_dir,'preprocessing')
        if not os.path.isdir(self.preprocess_dir):
            os.mkdir(self.preprocess_dir)
        
        if not os.path.isdir(os.path.join(self.preprocess_dir,'task-colors')):
            os.mkdir(os.path.join(self.preprocess_dir,'task-colors'))
        if not os.path.isdir(os.path.join(self.preprocess_dir,'task-letters')):
            os.mkdir(os.path.join(self.preprocess_dir,'task-letters'))
        if not os.path.isdir(os.path.join(self.preprocess_dir,'task-rsa')):
            os.mkdir(os.path.join(self.preprocess_dir,'task-rsa'))
    
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
        print('converting {} convert2bids_letter-color_main_{}.json'.format(self.mri_subject,self.subject))
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
        
        ###################
        # RENAME localizer runs: 'task-colorsloc' to 'task-colors'
        ###################
        dir_path = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')
        # loops through functional files
        for F in os.listdir(dir_path):
            # remove run from localizer event file names
            if ('task-colorsloc' in F):
                F_new = F.replace('task-colorsloc', 'task-colors') # old, new
                os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                print('old={} , new={}'.format(F,F_new))
            elif ('task-lettersloc' in F):
                F_new = F.replace('task-lettersloc', 'task-letters') # old, new
                os.rename(os.path.join(dir_path,F),os.path.join(dir_path,F_new)) # old,new
                print('old={} , new={}'.format(F,F_new))
                    
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
            behav_sess = 'sess-1' #old 
        else:
            behav_sess = 'sess-2'
        
        self.rename_logfiles(behav_sess)            # rename session and runs
        ## WARNING!! don't remove timestamps twice on same participant
        self.remove_timestamps_logfiles(behav_sess) # removes the trailing timestamps from the psychopy output
        self.copy_logfiles(behav_sess)              # copy logfiles to bids_raw folder with mri data
        print('success: housekeeping')
    
    def raw_copy(self,):
        """ Copy from bids_raw directory into derivaties folder for further processing/analysis.
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
    
    def preprocess_fsf(self,task):
        # Creates the FSF files for each subject's first level analysis        
        
        template_filename = os.path.join(self.template_dir,'preprocessing_template.fsf')
    
        markers = [
            '[$OUTPUT_PATH]', 
            '[$NR_TRS]',        # number of volumes
            '[$NR_VOXELS]',     # total number of voxels
            '[$FWHM]',
            '[$EPI_EECHO]',     # EPI DWELL TIME, EFFECTIVE ECHO SPACING
            '[$EPI_TE]',        # EPI echo time
            '[$INPUT_FILENAME]', # BOLD data
            '[$FMAP]',
            '[$FMAP_MAG_BRAIN]',
            '[$T1_BRAIN]',
            '[$MNI_BRAIN]'
        ]
    
        BOLD = os.path.join(self.deriv_dir,self.subject,self.session,'func','{}_{}_task-{}_bold.nii.gz'.format(self.subject,self.session,task))
        # calculate size of input data
        nii = nib.load(BOLD).get_data() # only do once 
        nr_trs = str(nii.shape[-1])
        nr_voxels = str(nii.size)
            
        FSF_filename = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_preprocessing_{}_{}.fsf'.format(task,self.subject,self.session) ) # save fsf
        # replacements
        output_path = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}'.format(task,self.subject,self.session)) 

        FMAP = os.path.join(self.deriv_dir,self.subject,self.session,'fmap','{}_{}_acq-fmap.nii.gz'.format(self.subject,self.session))
        FMAP_MAG_BRAIN = os.path.join(self.deriv_dir,self.subject,self.session,'fmap','{}_{}_run-01_fmap_brain.nii.gz'.format(self.subject,self.session))
        T1_BRAIN = os.path.join(self.deriv_dir,self.subject,self.session,'anat','{}_{}_T1w_brain.nii.gz'.format(self.subject,self.session))
        MNI_BRAIN  = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz')
        
        if task == 'rsa':
            FWHM = '0' # turn smoothing off for mulitvariate analyses
        else:
            FWHM = self.FWHM
            
        replacements = [ # needs to match order of 'markers'
            output_path,
            nr_trs,
            nr_voxels,
            FWHM,
            self.EPI_EECHO, # dwell time is effective echo spacing (EPI data not field map!!)
            self.EPI_TE,
            BOLD,
            FMAP,
            FMAP_MAG_BRAIN,
            T1_BRAIN,
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
        print('success: {}_fsf_preprocessing'.format(self.subject))
        
    def transform_2_mni(self, task, linear=1):
        """ Use the registration from the preprocessing FEAT output to transform the filtered_func_data.nii.gz to MNI space (non-linear)
        FLIRT and FNIRT registrations have already been calculated based on example_func.
        Here, we are just applying the existing FNIRT warps to the preprocessed data
        See here: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/FAQ#How_do_I_transform_a_mask_with_FLIRT_from_one_space_to_another.3F
        """
        
        # TIME SERIES TO BE TRANSFORMED
        EPI = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,self.session),'filtered_func_data.nii.gz')
        # EPI time series output non-linear NIFTI
        EPI_MNI = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,self.session),'filtered_func_data_mni.nii.gz')
        # standard space
        MNI = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz') # nifti
        
        if linear:
            # Apply FLIRT matrix (linear)
            example_func2standard = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,self.session),'reg','example_func2standard') 
            cmd = 'flirt -in {} -ref {} -applyxfm -init {}.mat -out {}'.format(EPI,MNI,example_func2standard,EPI_MNI)
            print(cmd)
            results = subprocess.call(cmd, shell=True, bufsize=0)
        else: 
            # Apply FNIRT warpfile (non-linear)
            # EPI time series (in MNI space) -> apply FNIRT warpfile based on T1 (in MNI space) -> warped to MNI space
            # commandline = 'applywarp -i input -o output -r reference -w warpfile/coefficients'
            example_func2standard_warp = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,self.session),'reg','example_func2standard_warp.nii.gz') 
            cmd = 'applywarp -i {} -o {} -r {} -w {}'.format(EPI,EPI_MNI,MNI,example_func2standard_warp)
            print(cmd)
            results = subprocess.call(cmd, shell=True, bufsize=0)
        print('success: transform_2_mni')
    
    def create_native_target(self, task):
        """ Copy the first session's mean_func to preprocessing folder and rename as native target
        """
        
        # use first session as registration target for both sessions
        native_target = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_reg_target.nii.gz'.format(task,self.subject))
        
        ###################
        # copy session 1's example func to task directory
        ###################
        src = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,'ses-01'),'reg','mean_func.nii.gz') # session1
        dst = reg_target
        sh.copyfile(src, dst)
        print('colors: src={} , dst={}'.format(src,dst))
        print('success: create_native_target')
        
    def transform_2_native_target(self, task):
        """ To  work in NATIVE space only, we need to create 1 target example_func for both sessions, FLIRT to the target, before concantenating EPIs.
        Use the first session's example_func as the native-space registration target for both sessions
        """
        # TIME SERIES TO BE TRANSFORMED
        EPI = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,self.session),'filtered_func_data.nii.gz')
        EPI_NATIVE = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,self.session),'filtered_func_data_native')
        native_target = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_native_target'.format(task,self.subject))
        
        # CREATE FLIRT transform (linear) on example_func
        example_func = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}.feat'.format(task,self.subject,self.session),'reg','example_func') 
        cmd = 'flirt -in {} -ref {}.nii.gz -out {}'.format(example_func,native_target,native_target) # save transforms in task folder preprocessing
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        # APPLY FLIRT transform (linear) to EPI
        cmd = 'flirt -in {} -ref {}.nii.gz -applyxfm -init {}.mat -out {}'.format(EPI,native_target,native_target,EPI_NATIVE) # what's the name of the matrix here?!
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)

        print('success: transform_2_native_target')