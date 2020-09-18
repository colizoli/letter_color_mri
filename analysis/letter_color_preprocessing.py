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

import os, subprocess, sys
import shutil as sh
import multiecho_combination as me
import nibabel as nib
from IPython import embed as shell

home_dir = '/nfs/colizoli/LetterColor'
analysis_dir = os.path.join(home_dir, 'analysis')
source_dir = os.path.join(home_dir, 'source') 
deriv_dir = os.path.join(home_dir,'derivatives')
    
def dicom2bids(subject, session, config):
    """ Converts MRI data from dicom to nifti format, then structures them according to the bids standard. 
    Takes about 15 minutes per session.
    Requirements:
        dcm2niix (apt get install dcm2niix)
        dcm2bids (python -m pip install dcm2bids) 
    """

    DICOM_DIR = os.path.join(source_dir,'dicom', subject, session)
    CONFIG_FILE  = os.path.join(analysis_dir, config)
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


def combine_echoes(subject, session, task, run, overwrite=False):
    """ TODO: Need to fix python version issue on lisa... either install as command line tool, or import packages to python3
    calling this tool from command line: https://github.com/dangom/multiecho
    run on TRIMMED echos if removing first N volumes for scanner stabilization
    JSON files need to match echo paths (so don't delete them)
    only runs in python 3.6
    """

    algorithm = 'te' # use echo time weighting, flat maps 

    output_name = os.path.join(deriv_dir,subj,session,'func','{}_task-{}_run-{}_bold.nii.gz'.format(subj,task,run))
    # Skip already combined images, in case it's rerun due to memory errors on some combines
    if not overwrite:
        if os.path.isfile(output_name):
            print('skipping subject {}, session {}, task {}, run {}'.format(subj,session,task,run))
            return
    
    this_run = os.path.join(deriv_dir,subj,session,'func','{}_task-{}_run-{}_echo-*_bold.nii.gz'.format(subj,task,run))
    
    # USE COMMAND LINE TOOL
    # cmd = 'mecombine \'{}\' --algorithm \'{}\' --outputname \'{}\''.format(this_run, algorithm, output_name) # command line
    # print(cmd)
    # results = subprocess.run(cmd, shell=True, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    # if results.returncode != 0:
    #     SubjectLog.subjects[subject].log('error_combine_echoes_{}_task-{}_run-{}'.format(session, task, run), results.stderr.decode("utf-8") + results.stdout.decode("utf-8"))
    #     return
    
    # USE AS PYTHON PACKAGE
    #echoes = me.load_me_data(this_run, None) # Called from inside me_combine() anyway
    combined, weights = me.me_combine(this_run, algorithm=algorithm) # combined returned as a nifti image with header and affine
    nib.save(combined, output_name) # data, fname

    # make sure there are no NaNs in the combined echos, zeros are generally better to work with
    cmd = 'fslmaths {} -nan {}'.format(output_name, output_name) # command line
    print(cmd)
    results = subprocess.call(cmd, shell=True, bufsize=0)
    
     

##### RUN PREPROCESSING #####
subjects = ['sub-503']
sessions = ['sess-01']
config = ['convert2bids_sub-503.json']


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
        dicom2bids(subj, session, config[s])
        source_copy(subj, session) # TODO don't forget to deface!

    for session in sessions:
        bet_brains_fmap(subj,session)   # B0 unwarping needs 'tight' brain extracted magnitude images of field map, better too small than too big!
        prepare_fmap(subj, session, echo_tdiff)  # prepares the field map image in radians/sec
        combine_echoes(subj, session, 'rsa', 2)  # for sub-503
    bet_brains_T1(subj) # brain extraction always check visually!