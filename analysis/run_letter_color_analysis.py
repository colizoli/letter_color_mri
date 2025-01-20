#!/usr/bin/env python
# encoding: utf-8

"""
letter_color_housekeeping.py
Created by O.Colizoli, June-2025
Last update: 09-01-2025
Python version 3.9

The following packages need to be installed because they are called from the command line, not imported:

fsl

"""

# TO DO:
# deface images T1s and fieldmaps before publication


# -----------------------
# Import functionality
# ----------------------- 
import os, subprocess, sys
import pandas as pd
import numpy as np
from IPython import embed as shell # for Oly's debugging only
# custom analysis scripts
import letter_color_housekeeping
import letter_color_first_level
# import letter_color_higher_level

# -----------------------
# Paths
# ----------------------- 
##########################
home_dir = os.path.dirname(os.getcwd()) # one level up from analysis folder
# home_dir = '/project/3018051.01/'
##########################
analysis_dir    = os.path.join(home_dir, 'analysis')        # scripts + configuration files for BIDS
# source_dir      = os.path.join(home_dir, 'raw')             # DICOMS DCCN project storage folder 'raw', directly imported from scanner (don't change!)
bids_dir        = os.path.join(home_dir, 'bids')        # NIFTI versions of DICOM files for processing
deriv_dir       = os.path.join(home_dir, 'derivatives')      # Processed data, fmriprep
mask_dir        = os.path.join(deriv_dir, 'masks')           # brain masks
template_dir    = os.path.join(deriv_dir, 'templates')    # fsl templates
timing_files_dir = os.path.join(deriv_dir, 'timing_files')   # custom 3 column format for 1st levels

# -----------------------
# Parameters
# -----------------------               
TR = 1.5 # seconds

# -----------------------
# Levels (switch ON/OFF)
# ----------------------- 
run_housekeeping    = False    # change file names, move files, etc.
run_first_level     = True     # concatenate runs, timing files, 1st level GLMs
run_higher_level    = False    # group-level analyses and statistics

# -----------------------
# Participants
# -----------------------
participants    = pd.read_csv(os.path.join(home_dir, 'participants_full_mri.csv'), dtype=str) # open in textmate, not excel!
subjects  = participants['subjects']


# -----------------------
# Run housekeeping class
# -----------------------

if run_housekeeping:
    for s,subject in enumerate(subjects): # go in order of mri        
        # initialize class
        housekeeping = letter_color_housekeeping.housekeeping_class(
            subject         = subject, # experiment subject number
            bids_dir        = bids_dir,
            deriv_dir       = deriv_dir,
            )            
        # housekeeping.rename_behav_logfiles()      # rename behavioral logfile names
        shell()
        
        # shell()
          
# -----------------------
# Run first-level class
# -----------------------   
if run_first_level:
    for s,subject in enumerate(subjects):
        # sessions are looped depending on task
        first_level = letter_color_first_level.first_level_class(
            subject         = subject,
            analysis_dir    = analysis_dir,
            bids_dir        = bids_dir,
            deriv_dir       = deriv_dir,
            mask_dir        = mask_dir,
            template_dir    = template_dir,
            timing_files_dir = timing_files_dir,
            TR              = TR, # repetition time in seconds
            )
        first_level.loc_match_bold()          # match loc1 and loc2 in nifti file to letters and colors in events files
        # first_level.loc_combine_epi('colors')             # concantenate both runs of localizer to perform 1 GLM
        # first_level.loc_combine_timing_files('colors')    # timing files for color localizer GLM
        # first_level.loc_combine_motion_regressors('colors')     # concatenate motion parameters from preprocessing, also outputs cols of 1s for each blocks' mean
        # first_level.loc_fsf('colors')                     # generates the first level FSF for the localizers
        # 
        # first_level.loc_combine_epi('letters')            # concantenate both runs of localizer to perform 1 GLM
        # first_level.loc_combine_timing_files('letters')   # timing files for color localizer GLM
        # first_level.loc_nuisance_regressors('letters')    # concatenate motion parameters from preprocessing, also outputs cols of 1s for each blocks' mean
        # first_level.loc_fsf('letters')                    # generates the first level FSF for the localizers
                
        # first_level.rsa_combine_epi()                     # concatenate EPI data for the 4 runs of the RSA task
        # first_level.rsa_combine_events()                  # concatenate events files for the 4 runs of the RSA task
        # first_level.rsa_dcm_split_nifti()                 # for spm, split the four-run concatenated nifti into single volume images and unzip
        # first_level.rsa_nuisance_regressors()             # volume-based physiological components, motion parameters, cosine (low-fres), and run means
        
        # first_level.rsa_timing_files_oddballs()           # create timing files for oddball stimuli
        # first_level.rsa_timing_files_letters()            # each letter in it's color and black: trained/untrained vs. color/black
        # first_level.rsa_timing_files_2x2()                # simple 2x2 design: trained/untrained vs. color/black
        
        # first_level.rsa_letters_fsf()                     # generates the first level FSF for the RSA design
        # first_level.rsa_2x2_fsf()                         # generates the first level FSF for the 2x2 design
        shell()
# -----------------------
# Run higher-level class
# -----------------------
if run_higher_level:
    ### PARTICIPANT LIST FOR BEHAVIOR ONLY
    participants = pd.read_csv(os.path.join(analysis_dir,'participants_full_behav.csv'), dtype=str) # open in textmate, not excel!
    higher_level = letter_color_higher_level.higher_level_class(
            subjects     = subjects_group,
            sessions     = sessions,
            analysis_dir = analysis_dir,
            deriv_dir    = deriv_dir,
            mask_dir     = mask_dir,
            template_dir = template_dir,
            TR           = TR, # repitition time in seconds
            participants = participants
            )
        
    higher_level.dataframe_trained_letters()                # all subjects letter conditions and colors 
    
    # higher_level.dataframe_qualia()                       # calculates the PA score from the questionnaire
    # higher_level.plot_qualia()                            # plots histograms and means of the PA scores
    
    # higher_level.dataframe_choose_pairs()                 # all group1 and group2's choices in a single dataframe
    # higher_level.plot_choose_pairs()                      # pie charts for frequency of color choices
    
    # higher_level.dataframe_subjects_iconic_memory()       # generates the dataframe for iconic memory higher level analysis
    # higher_level.dataframe_anova_iconic_memory()          # generates the dataframe for the ANOVA in JASP
    # higher_level.plot_anova_iconic_memory()                 # plot ANOVA 2x2 (collapsed over groups)

    # higher_level.dataframe_subjects_cieluv()              # converts from hexcodes to cieluv all subjects letters
    # higher_level.dataframe_subjects_consistency()         # generates the dataframe for the consistency score all subjects all letters
    # higher_level.dataframe_anova_consistency()            # generates the dataframe for the ANOVA in JASP
    # higher_level.plot_anova_consistency()                 # plot ANOVA 2x2 (collapsed over groups)
    
    # higher_level.correlation_consistency_iconic_memory()  # test correlation in the letter conditions between consistency scores and iconic memory performance
    
    # higher_level.localizers_randomise_input('letters')    # concatenates all subjects' cope1 in 4th dimension (localizers)
    # higher_level.localizers_randomise_input('colors')     # concatenates all subjects' cope1 in 4th dimension (localizers)
    
    # higher_level.rsa_letters_ev_conditions()      # outputs a DF with the letter and color conditions (general)
    # higher_level.rsa_letters_conditions()         # concatenates all subjects events files for letter-color conditions
    # higher_level.rsa_letters_combine_events()     # concatenates all subjects events files trial-wise
