#!/usr/bin/env python
# encoding: utf-8

"""
run_letter_color_analysis.py
Created by O.Colizoli, 2021
Last update: 21-2-2022
Python version 3.6

The following packages need to be installed because they are called from the command line, but not imported:
dcm2niix (conda install -c conda-forge dcm2niix)
dcm2bids (python -m pip install dcm2bids)
fsl

"""

# TO DO:
# Trim ends of RSA runs - sometimes scanner stops early, sometimes late
# sess-1 to ses-mri01 ?
# consider individual configs for subjects at bids level? Can specifiy which run is which localizer, and ignore bad scans (using series numbers?)
# cut physiological recordings into EPI time series
# implement RETROICOR

# PUBLICATION
# deface images T1s and fieldmaps before publication


# -----------------------
# Import functionality
# ----------------------- 
import os, subprocess, sys
import pandas as pd
from IPython import embed as shell # for Oly's debugging only
# custom analysis scripts
import letter_color_preprocessing
import letter_color_first_level

# -----------------------
# Paths
# ----------------------- 
##########################
# home_dir = os.path.dirname(os.getcwd()) # one level up from analysis folder
home_dir = '/project/3018051.01/'
##########################
analysis_dir = os.path.join(home_dir, 'analysis')       # scripts + configuration files for BIDS
source_dir = os.path.join(home_dir, 'raw')              # DICOMS DCCN project storage folder 'raw', directly imported from scanner (don't change!)
raw_dir = os.path.join(home_dir, 'bids_raw')            # NIFTI versions of DICOM files for processing
deriv_dir = os.path.join(home_dir,'derivatives')            # Processed data
mask_dir = os.path.join(deriv_dir,'masks')                  # brain masks
template_dir = os.path.join(analysis_dir,'templates')          # fsl templates
timing_files_dir = os.path.join(deriv_dir,'timing_files')   # custom 3 column format for 1st levels

# -----------------------
# Levels (switch ON/OFF)
# ----------------------- 
run_preprocessing = False    # motion correction, unwarping, registration, filtering, retroicor
run_first_level = True     # concatenate runs, timing files, 1st level GLMs
run_higher_level = False    # group-level analyses and statistics

# -----------------------
# Participants
# -----------------------
participants = pd.read_csv(os.path.join(analysis_dir,'participants.csv'), dtype=str) # open in textmate, not excel!
mri_subjects = participants['mri_subjects']
subjects_group = participants['subjects']
sessions = ['ses-01','ses-02']


# -----------------------
# Preprocessing class
# -----------------------
if run_preprocessing:
    for s,subject in enumerate(mri_subjects): # go in order of mri
        for session in sessions:
            # initialize class
            preprocess = letter_color_preprocessing.preprocess_class(
                subject = subjects_group[s], # experiment subject number
                mri_subject = subject,  # mri subject number
                session = session,
                analysis_dir = analysis_dir,
                source_dir = source_dir,
                raw_dir = raw_dir,
                deriv_dir = deriv_dir,
                mask_dir = mask_dir,
                template_dir = template_dir,
                timing_files_dir = timing_files_dir,
                )            
            preprocess.dicom2bids()       # convert DICOMS from scanner to nifti in bids format
            # preprocess.housekeeping()       # copies event files, rename file names to be bids compliant and same b/t mri & behavior
            # preprocess.raw_copy()         # copy from bids_raw directory into derivaties to prevent overwriting
            # preprocess.bet_brains_T1()    # brain extraction T1 always check visually!
            # preprocess.bet_brains_fmap()  # B0 unwarping needs 'tight' brain extracted magnitude images of field map, better too small than too big!
            # preprocess.prepare_fmap()     # prepares the field map image in radians/sec
            ### To-do!!
            # preprocess.preprocess_fsf()   # generate FSF file for preprocessing in FEAT
            # RETROICOR
            # preprocess.smooth_epi(fwhm=FWHM)       # smooths EPI data (localizers only)
        
# -----------------------
# First-level class
# -----------------------   
TR = 1.5 # in seconds
FWHM = 3 # smoothing kernel in mm (localizers)
# Values for unwarping in FEAT, based on BOLD data
EPI_TE = 0.0396*1000 # in ms
EPI_EECHO = 0.000580009*1000 # in ms

if run_first_level:
    for s,subject in enumerate(subjects_group):
        for session in sessions:
            first_level = letter_color_first_level.first_level_class(
                subject = subject,
                session = session,
                analysis_dir = analysis_dir,
                deriv_dir = deriv_dir,
                mask_dir = mask_dir,
                template_dir = template_dir,
                timing_files_dir = timing_files_dir,
                TR = TR, # repitition time in seconds
                FWHM = FWHM,
                EPI_TE = EPI_TE,
                EPI_EECHO = EPI_EECHO
                )
            # first_level.loc_timing_files('colors') # color localizer
            first_level.loc_fsf('colors')
        
        
            # first_level.concatenate_rsa_epi_data()           # concatenate EPI data for the 4 runs of the RSA task
            # first_level.rsa_timing_files('rsa', subject, session, rsa_runs[s]) # rsa, run as localizer color vs. black for pilot
            # # first_level.concatenate_timing_files()       # add time to block 2's onsets
            # # first_level.create_nuisance_regressors()     # concatenate motion parameters from preprocessing, also outputs cols of 1s for each blocks' mean
            # first_level.blocked_fsf()                     # create subject-specific FSF files for feat
    
    
# -----------------------
# Higher-level class
# -----------------------
    
# if higher_level:
    
    
