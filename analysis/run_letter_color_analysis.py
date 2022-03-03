#!/usr/bin/env python
# encoding: utf-8

"""
run_letter_color_analysis.py
Created by O.Colizoli, 2021
Last update: 21-2-2022
Python version 3.6

The following packages need to be installed because they are called from the command line, not imported:
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
run_preprocessing = True    # motion correction, unwarping, registration, filtering, retroicor
run_first_level = False     # concatenate runs, timing files, 1st level GLMs
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
# Values for unwarping in FEAT, based on BOLD data
EPI_TE = 0.0396*1000            # echo time EPI in ms
EPI_EECHO = 0.000580009*1000    # effective echo spacing EPI in ms
TR = 1.5    # repetition time in seconds
FWHM = 3    # smoothing kernel in mm (localizers only)

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
                EPI_TE = EPI_TE,
                EPI_EECHO = EPI_EECHO,
                FWHM = FWHM
                )            
            ## preprocess.dicom2bids()       # DONE  till sub-226!!! convert DICOMS from scanner to nifti in bids format
            ## preprocess.housekeeping()     # DONE  till sub-226!!! till sub-226 copies event files, rename file names to be bids compliant and same b/t mri & behavior
            ## preprocess.raw_copy()         # DONE  till sub-226!!! copy from bids_raw directory into derivaties to prevent overwriting
            # preprocess.bet_brains_T1()    # brain extraction T1 always check visually!
            # preprocess.bet_brains_fmap()  # B0 unwarping needs 'tight' brain extracted magnitude images of field map, better too small than too big!
            # preprocess.prepare_fmap()     # prepares the field map image in radians/sec            
            
            try:
                # preprocess.preprocess_fsf('letters')          # generate FSF file for preprocessing in FEAT (run from command line - batch)
                # preprocess.preprocess_fsf('colors')           # generate FSF file for preprocessing in FEAT (run from command line - batch)
                # preprocess.preprocess_fsf('rsa','_run-01')    # generate FSF file for preprocessing in FEAT (run from command line - batch)
                # preprocess.preprocess_fsf('rsa','_run-02')    # generate FSF file for preprocessing in FEAT (run from command line - batch)
                # preprocess.preprocess_fsf('rsa','_run-03')    # generate FSF file for preprocessing in FEAT (run from command line - batch)
                # preprocess.preprocess_fsf('rsa','_run-04')    # generate FSF file for preprocessing in FEAT (run from command line - batch)
                
                preprocess.transform_2_mni('rsa','_run-01')     # transforms the preprocessed time series to MNI space
                preprocess.transform_2_mni('rsa','_run-02')     # transforms the preprocessed time series to MNI space
                preprocess.transform_2_mni('rsa','_run-03')     # transforms the preprocessed time series to MNI space
                preprocess.transform_2_mni('rsa','_run-04')     # transforms the preprocessed time series to MNI space
                preprocess.transform_2_mni('letters')           # transforms the preprocessed time series to MNI space
                preprocess.transform_2_mni('colors')            # transforms the preprocessed time series to MNI space

                preprocess.transform_2_native_target('rsa','_run-01')   # register to session 1 native space
                preprocess.transform_2_native_target('rsa','_run-02')   # register to session 1 native space
                preprocess.transform_2_native_target('rsa','_run-03')   # register to session 1 native space
                preprocess.transform_2_native_target('rsa','_run-04')   # register to session 1 native space
                preprocess.transform_2_native_target('letters')         # register to session 1 native space RSA task
                preprocess.transform_2_native_target('colors')          # register to session 1 native space
            except:
                pass
            

            ### To-do!!
            # RETROICOR
            
        
# -----------------------
# First-level class
# -----------------------   

if run_first_level:
    for s,subject in enumerate(subjects_group):
        # sessions are looped depending on task
        first_level = letter_color_first_level.first_level_class(
            subject = subject,
            analysis_dir = analysis_dir,
            deriv_dir = deriv_dir,
            mask_dir = mask_dir,
            template_dir = template_dir,
            timing_files_dir = timing_files_dir,
            TR = TR, # repitition time in seconds
            )
        # first_level.loc_combine_epi('colors')    # concantenate both runs of localizer to perform 1 GLM
        # first_level.loc_combine_timing_files('colors')    # timing files for color localizer GLM
        # first_level.loc_nuisance_regressors('colors')     # concatenate motion parameters from preprocessing, also outputs cols of 1s for each blocks' mean
        # first_level.loc_fsf('colors',run_cmd=0)           # generates the first level FSF for the localizers
        
        # first_level.loc_combine_epi('letters')    # concantenate both runs of localizer to perform 1 GLM
        # first_level.loc_combine_timing_files('letters')   # timing files for color localizer GLM
        # first_level.loc_nuisance_regressors('letters')    # concatenate motion parameters from preprocessing, also outputs cols of 1s for each blocks' mean
        # first_level.loc_fsf('letters',run_cmd=0)          # generates the first level FSF for the localizers
        
        # first_level.rsa_combine_epi()                     # concatenate EPI data for the 4 runs of the RSA task
        # first_level.rsa_combine_events()                  # concatenate events files for the 4 runs of the RSA task
        # first_level.rsa_nuisance_regressors()             # motion parameters, run means, and oddball trials as nuisance
        
        # first_level.rsa_timing_files_2x2()                # simple 2x2 design: trained/untrained vs. color/black
        # first_level.rsa_2x2_fsf(run_cmd=0)                # generates the first level FSF for the 2x2 design
        
        # TO DO 
        # first_level.rsa_timing_files_letters()                # simple 2x2 design: trained/untrained vs. color/black
        # first_level.rsa_letters_fsf(run_cmd=1)                # generates the first level FSF for the RSA design
        
    
# -----------------------
# Higher-level class
# -----------------------
    
# if higher_level:
    
    
