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
import numpy as np
from IPython import embed as shell # for Oly's debugging only

class preprocess_class(object):
    def __init__(self, subject, mri_subject, session, analysis_dir, source_dir, raw_dir, deriv_dir, mask_dir, template_dir, FWHM, t1_path):        
        self.subject        = 'sub-'+str(subject)
        self.mri_subject    = 'sub-'+str(mri_subject)
        self.session        = str(session)
        self.analysis_dir   = str(analysis_dir)
        self.source_dir     = str(source_dir)   # dicoms
        self.raw_dir        = str(raw_dir)      # raw niftis
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.FWHM           = str(FWHM)
        self.t1_path        = str(t1_path) # with _brain extension
        
        # path to registration target for both sessions (output)
        self.native_target = os.path.join(self.preprocess_dir,'{}_native_target.nii.gz'.format(self.subject))
            
        if not os.path.isdir(self.raw_dir):
            os.mkdir(self.raw_dir)
        if not os.path.isdir(self.deriv_dir):
            os.mkdir(self.deriv_dir)
        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
        
        self.preprocess_dir = os.path.join(self.deriv_dir,'preprocessing')
        if not os.path.isdir(self.preprocess_dir):
            os.mkdir(self.preprocess_dir)
        
        if not os.path.isdir(os.path.join(self.preprocess_dir,'task-colors')):
            os.mkdir(os.path.join(self.preprocess_dir,'task-colors'))
        if not os.path.isdir(os.path.join(self.preprocess_dir,'task-letters')):
            os.mkdir(os.path.join(self.preprocess_dir,'task-letters'))
        if not os.path.isdir(os.path.join(self.preprocess_dir,'task-rsa')):
            os.mkdir(os.path.join(self.preprocess_dir,'task-rsa'))
        
        # write unix commands to job to run in parallel
        self.preprocessing_job_path = os.path.join(self.analysis_dir,'jobs','job_preprocessing_{}.txt'.format(self.subject))
        if self.session == 'ses-01':
            self.preprocessing_job = open(self.preprocessing_job_path, "w")
            self.preprocessing_job.write("#!/bin/bash\n")
            self.preprocessing_job.close()
    
    def dicom2bids(self, ):
        """ Convert MRI data from dicom to nifti format, then structures them according to the bids standard. 
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
        """ Renames the log files to match MRI files: e.g., sess-2 -> ses-02
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
        """ removes the trailing timestamps from the psychopy output
         e.g. last 16 characters _00200312-125726
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
        for f in os.listdir(os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')):
            src = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func', f)
            dst = os.path.join(self.raw_dir, 'nifti', self.subject, self.session, 'func', f)
            sh.copy(src, dst)
        print('func: src={} , dst={}'.format(src,dst))
        
        print('success: copy_logfiles')
    
    def copy_physio(self, behav_sess):
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
        for f in os.listdir(os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func')):
            src = os.path.join(self.raw_dir, 'logfiles', self.subject, behav_sess, 'func', f)
            dst = os.path.join(self.raw_dir, 'nifti', self.subject, self.session, 'func', f)
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
        
        # self.rename_logfiles(behav_sess)            # rename session and runs
        # ## WARNING!! don't remove timestamps twice on same participant
        # self.remove_timestamps_logfiles(behav_sess) # removes the trailing timestamps from the psychopy output
        # self.copy_logfiles(behav_sess)              # copy logfiles to bids_raw folder with mri data
        print('success: housekeeping')
    
    
    def raw_copy(self,):
        """ Copy from bids_raw directory into derivatives folder for further processing/analysis.
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
    
    
    def create_native_target(self, task='rsa', session='ses-01', bold_run='run-03'):
        """Create the egistration target in native space for ALL tasks. 
        
        Args:
            task (str): rsa or localizer task to use for target
            session (str): ses-01 or ses-02 to use for target
            bold_run (str): the bold run to use for target
        
        Notes:
        ------
            The target is the first session, 3rd run (RSA3), motion corrected and then mean image.
            Path is defined at class level: self.native_target
        """        
        # path to raw bold file as input to create registration target
        mri_in = os.path.join(self.deriv_dir,self.subject,self.session,'func','{}_{}_task-{}_{}_bold.nii.gz'.format(self.subject, self.session, task, bold_run))
        mri_out = os.path.join(self.preprocess_dir, '{}_{}_task-{}_{}_bold_mcf.nii.gz'.format(self.subject, self.session, task, bold_run))
        
        ###################
        # motion correction (default target middle volume)
        ###################
        cmd = 'mcflirt -in {} -stages 4 -sinc_final -o {}'.format(mri_in, mri_out)
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        
        ###################
        # mean image
        ###################
        cmd = 'fslmaths {} -Tmean {}'.format(mri_out, self.native_target)
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        
        ###################
        # brain extraction
        ###################
        cmd = 'bet {} {}'.format(self.native_target, '_brain')
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        

        # open preprocessing job and write command as new line
        # cmd = 'cp {} {}'.format(src,dst)
        # self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
        # self.preprocessing_job.write(cmd)   # command
        # self.preprocessing_job.write("\n\n")  # new line
        # self.preprocessing_job.close()
        print('success: create_native_target')
    
    
    def native_target_2_mni(self, ):
        """Compute the registration from the subject-specific native target to MNI space via the T1.
        
        Notes:
        ------
        FLIRT to T1 then FNIRT to MNI.
        Save transformation matrices to apply later.
        """
        # path to registration target for both sessions (output)
        # split self.native target and add _mni to end
        mri_out = os.path.join(self.preprocess_dir,'{}_native_target_mni.nii.gz'.format(self.subject))
        
        example_func = self.native_target
        
        ###################
        # FLIRT to T1
        ###################
        # /usr/local/fsl/bin/epi_reg --epi=example_func --t1=highres_head --t1brain=highres --out=example_func2highres
        #
        # /usr/local/fsl/bin/convert_xfm -inverse -omat highres2example_func.mat example_func2highres.mat
        #
        # /usr/local/fsl/bin/flirt -in highres -ref standard -out highres2standard -omat highres2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear
        #
        # /usr/local/fsl/bin/fnirt --iout=highres2standard_head --in=highres_head --aff=highres2standard.mat --cout=highres2standard_warp --iout=highres2standard --jout=highres2highres_jac --config=T1_2_MNI152_2mm --ref=standard_head --refmask=standard_mask --warpres=10,10,10
        #
        # /usr/local/fsl/bin/applywarp -i highres -r standard -o highres2standard -w highres2standard_warp
        #
        # /usr/local/fsl/bin/convert_xfm -inverse -omat standard2highres.mat highres2standard.mat
        #
        # /usr/local/fsl/bin/convert_xfm -omat example_func2standard.mat -concat highres2standard.mat example_func2highres.mat
        #
        # /usr/local/fsl/bin/convertwarp --ref=standard --premat=example_func2highres.mat --warp1=highres2standard_warp --out=example_func2standard_warp
        #
        # /usr/local/fsl/bin/applywarp --ref=standard --in=example_func --out=example_func2standard --warp=example_func2standard_warp
        #
        # /usr/local/fsl/bin/convert_xfm -inverse -omat standard2example_func.mat example_func2standard.mat
        
        
        
        cmd = 'flirt -in {} -o {}'.format(self.native_target, mri_out)
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        
        ###################
        # FNIRT to MNI
        ###################
        cmd = 'fslmaths {} -Tmean {}'.format(mri_out, native_target)
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        
        
        # open preprocessing job and write command as new line
        # self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
        # self.preprocessing_job.write(cmd)   # command
        # self.preprocessing_job.write("\n\n")  # new line
        # self.preprocessing_job.close()
        print('success: native_target_2_mni')
    
    
    def motion_correction(self, ):
        """Run motion correction with subject-specific native target as reference volume.
        
        Notes:
        ------
        Call MCFLIRT from command line (FEAT has no option for overriding default reference volume.)
        Save motion matrices to add to first level GLM.
        """
        # go into subject's session's func folder and motion correct all nii.gz files
        for task in ['letters','colors','rsa_run-01','rsa_run-02','rsa_run-03','rsa_run-04']:
            # raw bold data 
            mri_in = os.path.join(self.deriv_dir,self.subject,self.session,'func','{}_{}_task-{}_bold.nii.gz'.format(self.subject,self.session,task))
            
            if os.path.exists(mri_in): ##### SKIP IF BOLD FILE DOES NOT EXIST #####
                # save in preprocessing folder of task
                mri_out = os.path.join(self.preprocess_dir,'task-{}'.format(task),'{}_{}_task-{}_bold_mcf'.format(self.subject,self.session,task))
        
                ###################
                # motion correction (default target middle volume)
                ###################
                cmd = 'mcflirt -in {} -r {} -stages 4 -sinc_final -o {} -stats -mats -plots'.format(mri_in, self.native_target, mri_out)
                print(cmd)
                results = subprocess.call(cmd, shell=True, bufsize=0)
                # # open preprocessing job and write command as new line
                # self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                # self.preprocessing_job.write(cmd)   # command
                # self.preprocessing_job.write("\n\n")  # new line
                # self.preprocessing_job.close()                    
        print('success: motion_correction')
        
        
    def preprocess_fsf(self, ):
        # Creates the FSF files for each subject's first level analysis        
        
        template_filename = os.path.join(self.template_dir,'preprocessing_template.fsf')
        
        markers = [
            '[$OUTPUT_PATH]', 
            '[$NR_TRS]',        # number of volumes
            '[$NR_VOXELS]',     # total number of voxels
            '[$FWHM]',
            '[$INPUT_FILENAME]', # BOLD data
            '[$T1_BRAIN]',
            '[$MNI_BRAIN]'
        ]
        
        # all nii.gz files in func folder need to be preprocessed in FEAT
        for task in ['letters','colors','rsa_run-01','rsa_run-02','rsa_run-03','rsa_run-04']:
            
            # This will be the mcflirted time series in the preprocessing folder
            mri_in = os.path.join(self.preprocess_dir,'task-{}'.format(task),'{}_{}_task-{}{}_bold_mcf.nii.gz'.format(self.subject,self.session,task,bold_run))

            if os.path.exists(mri_in):  ##### SKIP IF BOLD FILE DOES NOT EXIST #####
                # calculate size of input data
                nii = nib.load(mri_in).get_data() # only do once 
                NR_TRS = str(nii.shape[-1])
                NR_VOXELS = str(nii.size)
        
                FSF_filename = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_preprocessing_{}_{}{}.fsf'.format(task,self.subject,self.session,bold_run) ) # save fsf
                # replacements
                output_path = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}{}'.format(task,self.subject,self.session,bold_run)) 

                t1_brain = self.t1_path
                mni_brain  = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz')
        
                if task == 'rsa':
                    FWHM = '0' # turn smoothing off for mulitvariate analyses
                else:
                    FWHM = self.FWHM
            
                replacements = [ # needs to match order of 'markers'
                    output_path,
                    NR_TRS,
                    NR_VOXELS,
                    FWHM,
                    mri_in,
                    t1_brain,
                    mni_brain
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
        
                # open preprocessing job and write command as new line
                cmd = 'feat {}'.format(FSF_filename)
                self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
                self.preprocessing_job.write(cmd)   # feat command
                self.preprocessing_job.write("\n\n")  # new line
                self.preprocessing_job.close()
                print('success: {}'.format(FSF_filename))
            else:
                print('cannot make FSF: missing {}'.format(BOLD))
        
    def transform_2_mni(self, task, bold_run='', linear=0):
        """ Use the registration from the preprocessing FEAT output to transform the filtered_func_data.nii.gz to MNI space (non-linear)
        FLIRT and FNIRT registrations have already been calculated based on example_func.
        Here, we are just applying the existing FNIRT warps to the preprocessed data
        See here: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/FAQ#How_do_I_transform_a_mask_with_FLIRT_from_one_space_to_another.3F
        """
        
        ##### SKIP IF BOLD INPUT FILE DOES NOT EXIST #####
        BOLD = os.path.join(self.deriv_dir,self.subject,self.session,'func','{}_{}_task-{}{}_bold.nii.gz'.format(self.subject,self.session,task,bold_run))
        if os.path.exists(BOLD):
            # TIME SERIES TO BE TRANSFORMED
            EPI = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}{}.feat'.format(task,self.subject,self.session,bold_run),'filtered_func_data.nii.gz')
            # EPI time series output non-linear NIFTI
            EPI_MNI = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}{}.feat'.format(task,self.subject,self.session,bold_run),'filtered_func_data_mni.nii.gz')
            # standard space
            MNI = os.path.join(self.mask_dir, 'MNI152_T1_2mm_brain.nii.gz') # nifti
        
            if linear:
                # Apply FLIRT matrix (linear)
                example_func2standard = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}{}.feat'.format(task,self.subject,self.session,bold_run),'reg','example_func2standard') 
                cmd = 'flirt -in {} -ref {} -applyxfm -init {}.mat -out {}'.format(EPI,MNI,example_func2standard,EPI_MNI)
                print(cmd)
                # results = subprocess.call(cmd, shell=True, bufsize=0)
            else: 
                # Apply FNIRT warpfile (non-linear)
                # EPI time series (in MNI space) -> apply FNIRT warpfile based on T1 (in MNI space) -> warped to MNI space
                # commandline = 'applywarp -i input -o output -r reference -w warpfile/coefficients'
                example_func2standard_warp = os.path.join(self.preprocess_dir,'task-{}'.format(task),'task-{}_{}_{}{}.feat'.format(task,self.subject,self.session,bold_run),'reg','example_func2standard_warp.nii.gz') 
                cmd = 'applywarp -i {} -o {} -r {} -w {}'.format(EPI,EPI_MNI,MNI,example_func2standard_warp)
                print(cmd)
                # results = subprocess.call(cmd, shell=True, bufsize=0)
            # open preprocessing job and write command as new line
            self.preprocessing_job = open(self.preprocessing_job_path, "a") # append is important, not write
            self.preprocessing_job.write(cmd)   # command
            self.preprocessing_job.write("\n\n")  # new line
            self.preprocessing_job.close()
        print('success: transform_2_mni')
