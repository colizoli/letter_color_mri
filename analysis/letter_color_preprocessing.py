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
    def __init__(self, subject, session, analysis_dir, source_dir,raw_dir, deriv_dir, mask_dir, template_dir, timing_files_dir):        
        self.subject        = 'sub-'+str(subject)
        self.session        = str(session)
        self.analysis_dir   = str(analysis_dir)
        self.source_dir     = str(source_dir)
        self.raw_dir        = str(raw_dir)
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.timing_files_dir = str(timing_files_dir)
        
        if not os.path.isdir(self.source_dir):
            os.mkdir(self.source_dir)
        if not os.path.isdir(self.deriv_dir):
            os.mkdir(self.deriv_dir)
        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
        if not os.path.isdir(self.timing_files_dir):
            os.mkdir(self.timing_files_dir)    

    
    def dicom2bids(self, config):
        """ Converts MRI data from dicom to nifti format, then structures them according to the bids standard. 
        Takes about 15 minutes per session.
        Requirements:
            dcm2niix (conda install -c conda-forge dcm2niix)
            dcm2bids (python -m pip install dcm2bids) 
        """

        DICOM_DIR = os.path.join(self.source_dir, self.subject, self.session) #/Users/olympiacolizoli/Aeneas/mountpoint6/raw
        CONFIG_FILE  = os.path.join(self.analysis_dir, 'config', config)
        OUTPUT_DIR = os.path.join(self.raw_dir, 'nifti', self.subject, self.session)
        # Create folder structure if not yet existing
        if not os.path.isdir(os.path.join(self.raw_dir, 'nifti')):
            os.mkdir(os.path.join(self.raw_dir, 'nifti'))
        if not os.path.isdir(os.path.join(self.raw_dir, 'nifti', self.subject)):
            os.mkdir(os.path.join(self.raw_dir, 'nifti', self.subject))
        if not os.path.isdir(os.path.join(self.raw_dir, 'nifti', self.subject, self.session)):
            os.mkdir(os.path.join(self.raw_dir, 'nifti', self.subject, self.session)) 
    
        cmd = 'dcm2bids -d {} -p {} -c {} -o {}'.format(DICOM_DIR, self.subject, CONFIG_FILE, OUTPUT_DIR)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        print('success: dicom2bids')


    def raw_copy(self,):
        """ copy from bids_raw directory into derivaties folder for further processing/analysis.
        All further analysis should be run within the derivatives folder.
        """

        # Create folder structure if not yet existing
        if not os.path.isdir(os.path.join(self.deriv_dir)):
            os.mkdir(os.path.join(self.deriv_dir))
        if not os.path.isdir(os.path.join(self.deriv_dir, self.subject)):
            os.mkdir(os.path.join(self.deriv_dir, self.subject))
        if os.path.isdir(os.path.join(self.deriv_dir, self.subject, self.session)):
            print('Target folder already exists. No files will be copied. Delete existing session from derivatives if intended.')
            return
        # Copy nifti
        src = os.path.join(self.raw_dir, 'nifti', self.subject, self.session, self.subject)
        dst = os.path.join(self.deriv_dir, self.subject, self.session)
        sh.copytree(src, dst)
        # Copy event files
        # src = os.path.join(self.raw_dir, 'events', self.subject, self.session)
        # dst = os.path.join(self.deriv_dir, self.subject, self.session, 'events')
        # sh.copytree(src, dst)
        print('success: source_copy')

    def bet_brains_T1(self, postfix='brain'):
        """ Runs brain extraction on all T1s.
        Always check visually in fsleyes.
        """
    
        inFile = os.path.join(self.deriv_dir, self.subject, 'sess-01', 'anat', '{}_T1w.nii.gz'.format(self.subject))
        outFile = os.path.join(self.deriv_dir, self.subject, 'sess-01', 'anat', '{}_T1w_{}.nii.gz'.format(self.subject, postfix))
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
            inFile = os.path.join(self.deriv_dir,  self.subject,  self.session, 'fmap', '{}_fmap_run-{}_fmap.nii.gz'.format(self.subject,run))
            outFile = os.path.join(self.deriv_dir,  self.subject,  self.session, 'fmap', '{}_fmap_run-{}_fmap_{}.nii.gz'.format(self.subject,run,postfix))
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
        'dwell_time' refers to effective echo spacing of EPI data
        echo_tdiff = 0.00246 # seconds, difference between echo times of field map magnitude images [0.0047,0.00716]
        """
    
        echo_tdiff = echo_tdiff * 1000 # needs to be in milliseconds
        phase_image = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_fmap_run-03_fmap.nii.gz'.format(self.subject))
        mag_image = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_fmap_run-01_fmap_brain.nii.gz'.format(self.subject))
        outFile = os.path.join(self.deriv_dir, self.subject, self.session, 'fmap', '{}_acq-fmap.nii.gz'.format(self.subject))
        # bet inFile outFile
        cmd = 'fsl_prepare_fieldmap {} {} {} {} {} [--nocheck]'.format('SIEMENS',phase_image,mag_image,outFile,echo_tdiff)
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
        print('success: prepare_fmap')
