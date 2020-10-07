#!/usr/bin/env python
# encoding: utf-8

"""
letter_color_preprocessing.py
Created by O.Colizoli, June-2019
Last update: 11-06-2019
Python version 3.4

The following packages need to be installed because they are called from the command line, but not imported:
Python3.6 with https://github.com/dangom/multiecho
dcm2niix
dcm2bids (installed through pip)
fsl

module load Python/3.6.3-foss-2017b
"""

# TO DO:
# Trim ends of RSA runs - sometimes scanner stops early, sometimes late
# sess-1 to ses-mri01 ?
# consider individual configs for subjects at bids level? Can specifiy which run is which localizer, and ignore bad scans (using series numbers?)

import os, subprocess, sys
import shutil as sh
# import multiecho_combination as me
import nibabel as nib
from IPython import embed as shell

##########################
home_dir = '/Users/olympiacolizoli/Aeneas/mountpoint6/'
##########################
analysis_dir = os.path.join(home_dir, 'analysis') # scripts + configuration files for BIDS
raw_dir = os.path.join(home_dir, 'raw') # DCCN project storage folder 'raw', directly imported from scanner (don't change!)
source_dir = os.path.join(home_dir, 'source') # NIFTI versions of DICOM files for processing
deriv_dir = os.path.join(home_dir,'derivatives') # Processed data
    
def dicom2bids(subject, session, config):
    """ Converts MRI data from dicom to nifti format, then structures them according to the bids standard. 
    Takes about 15 minutes per session.
    Requirements:
        dcm2niix (conda install -c conda-forge dcm2niix)
        dcm2bids (python -m pip install dcm2bids) 
    """

    DICOM_DIR = os.path.join(raw_dir, subject, session) #/Users/olympiacolizoli/Aeneas/mountpoint6/raw
    CONFIG_FILE  = os.path.join(analysis_dir, 'config', config)
    OUTPUT_DIR = os.path.join(source_dir, 'nifti', subject, session)
    # Create folder structure if not yet existing
    if not os.path.isdir(os.path.join(source_dir, 'nifti')):
        os.mkdir(os.path.join(source_dir, 'nifti'))
    if not os.path.isdir(os.path.join(source_dir, 'nifti', subject)):
        os.mkdir(os.path.join(source_dir, 'nifti', subject))
    if not os.path.isdir(os.path.join(source_dir, 'nifti', subject, session)):
        os.mkdir(os.path.join(source_dir, 'nifti', subject, session)) 
    
    cmd = 'dcm2bids -d {} -p {} -c {} -o {}'.format(DICOM_DIR, subject, CONFIG_FILE, OUTPUT_DIR)
    results = subprocess.call(cmd, shell=True, bufsize=0)
    print('success: dicom2bids')


def source_copy(subject, session):
    """ Add a copy to the derivatives folder for further processing/analysis.
    All further analysis should be run within the derivatives folder.
    """

    # Create folder structure if not yet existing
    if not os.path.isdir(os.path.join(deriv_dir)):
        os.mkdir(os.path.join(deriv_dir))
    if not os.path.isdir(os.path.join(deriv_dir, subject)):
        os.mkdir(os.path.join(deriv_dir, subject))
    if os.path.isdir(os.path.join(deriv_dir, subject, session)):
        print('Target folder already exists. No files will be copied. Delete existing session from derivatives if intended.')
        return
    # Copy nifti
    src = os.path.join(source_dir, 'nifti', subject, session, subject)
    dst = os.path.join(deriv_dir, subject, session)
    sh.copytree(src, dst)
    # Copy event files
    src = os.path.join(source_dir, 'events', subject, session)
    dst = os.path.join(deriv_dir, subject, session, 'events')
    sh.copytree(src, dst)
    print('success: source_copy')

def bet_brains_T1(subject, postfix='brain'):
    """ Runs brain extraction on all T1s.
    Always check visually in fsleyes.
    """
    
    inFile = os.path.join(deriv_dir, subject, 'sess-01', 'anat', '{}_T1w.nii.gz'.format(subj))
    outFile = os.path.join(deriv_dir, subject, 'sess-01', 'anat', '{}_T1w_{}.nii.gz'.format(subj, postfix))
    # bet inFile outFile
    cmd = 'bet {} {}'.format(inFile, outFile)
    print(cmd)
    results = subprocess.call(cmd, shell=True, bufsize=0)

    # reorient to mni space if necessary (sometimes turned around)
    # cmd = 'fslreorient2std {} {}'.format(outFile,outFile)
    # proc = subprocess.call( cmd, shell=True, bufsize=0,) # COMMAND LINE
    print('success: bet_brains_T1')
    
def bet_brains_fmap(subject, session, postfix='brain'):
    # Runs brain extraction on magnitude images of field maps
    # This needs to be 'tight', meaning it is better to exclude brain voxels than include noisy non-brain voxels
    # Important! Always check visually in fsleyes

    for run in ['01','02','03']: # first 2 are magnitude, 3rd is phase
        inFile = os.path.join(deriv_dir, subject, session, 'fmap', '{}_fmap_run-{}_fmap.nii.gz'.format(subj,run))
        outFile = os.path.join(deriv_dir, subject, session, 'fmap', '{}_fmap_run-{}_fmap_{}.nii.gz'.format(subj,run,postfix))
        # bet inFile outFile
        cmd = 'bet {} {} -f 0.6'.format(inFile, outFile) 
        print(cmd)
        results = subprocess.call(cmd, shell=True, bufsize=0)
    
        # reorient to mni space if necessary (sometimes turned around)
        # cmd = 'fslreorient2std {} {}'.format(outFile,outFile)
        # proc = subprocess.call( cmd, shell=True, bufsize=0,) # COMMAND LINE
    print('success: bet_brains_fmap')
    
def prepare_fmap(subject, session, echo_tdiff):
    """ Need to prepare field map before B0 unwarping
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide#SIEMENS_data
    Run after brain extraction on magnitude image
    Note too important which magnitude image is chosen
    Echo time difference need to be in ms (echotimes in JSON)
    fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE (in ms)> [--nocheck]
    """
    
    echo_tdiff = echo_tdiff * 1000 # needs to be in milliseconds
    phase_image = os.path.join(deriv_dir, subject, session, 'fmap', '{}_fmap_run-03_fmap.nii.gz'.format(subj))
    mag_image = os.path.join(deriv_dir, subject, session, 'fmap', '{}_fmap_run-01_fmap_brain.nii.gz'.format(subj))
    outFile = os.path.join(deriv_dir, subject, session, 'fmap', '{}_acq-fmap.nii.gz'.format(subj))
    # bet inFile outFile
    cmd = 'fsl_prepare_fieldmap {} {} {} {} {} [--nocheck]'.format('SIEMENS',phase_image,mag_image,outFile,echo_tdiff)
    print(cmd)
    results = subprocess.call(cmd, shell=True, bufsize=0)
    print('success: prepare_fmap')



##### RUN PREPROCESSING #####
# subjects = ['sub-101','sub-201','sub-102','sub-202','sub-103','sub-203','sub-104','sub-204','sub-105','sub-205','sub-106','sub-206','sub-107']
subjects = ['sub-001']
sessions = ['ses-mri01']
# sessions = ['ses-mri01','ses-mri02']
config = [ # order of localizers: LC = letters, colors; CL = colors, letters
    'convert2bids_letter-color_main_LC.json', # sub-001  
]


# TR4 = 1.5
# TR6 = 1.0
# rmTRs = 5 # TRs to be removed for stabilization of magnetic field
#### B0 unwarping ####
# current issue: FEAT asks for echo time to run unwarping, but what is the echo time for the combined images?
# however, when unwarping run from FUGUE, it only asks for dwell time and magnitude image.
# what is FEAT doing differently from FUGUE? (combining with motion correct and BBR, but why is echo time needed?)
# 'dwell_time' refers to effective echo spacing of EPI data
echo_tdiff = 0.00246 # seconds, difference between echo times of field map magnitude images [0.0047,0.00716]

########################

## LOOP SUBJECTS    
for s,subj in enumerate(subjects):
    # convert from dicom to difti, copy into derivaties to prevent overwriting
    for session in sessions:
        dicom2bids(subj, session, configs[s])
    #     source_copy(subj, session) # TODO don't forget to deface!
    #
    # for session in sessions:
    #     bet_brains_fmap(subj,session)   # B0 unwarping needs 'tight' brain extracted magnitude images of field map, better too small than too big!
    #     prepare_fmap(subj, session, echo_tdiff)  # prepares the field map image in radians/sec
    # bet_brains_T1(subj) # brain extraction always check visually!