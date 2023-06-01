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
# ?? Trim ends of RSA runs - sometimes scanner stops early, sometimes late

# PUBLICATION
# deface images T1s and fieldmaps before publication


# -----------------------
# Import functionality
# ----------------------- 
import os, subprocess, sys
import pandas as pd
import numpy as np
from IPython import embed as shell # for Oly's debugging only
# custom analysis scripts
import letter_color_preprocessing
import letter_color_first_level
import letter_color_higher_level

# -----------------------
# Paths
# ----------------------- 
##########################
# home_dir = os.path.dirname(os.getcwd()) # one level up from analysis folder
home_dir = '/project/3018051.01/'
##########################
analysis_dir    = os.path.join(home_dir, 'analysis')        # scripts + configuration files for BIDS
source_dir      = os.path.join(home_dir, 'raw')             # DICOMS DCCN project storage folder 'raw', directly imported from scanner (don't change!)
raw_dir         = os.path.join(home_dir, 'bids_raw')        # NIFTI versions of DICOM files for processing
deriv_dir       = os.path.join(home_dir,'derivatives')      # Processed data
mask_dir        = os.path.join(deriv_dir,'masks')           # brain masks
template_dir    = os.path.join(analysis_dir,'templates')    # fsl templates
timing_files_dir = os.path.join(deriv_dir,'timing_files')   # custom 3 column format for 1st levels

# -----------------------
# Levels (switch ON/OFF)
# ----------------------- 
run_preprocessing   = True    # motion correction, unwarping, registration, filtering, retroicor
run_first_level     = False     # concatenate runs, timing files, 1st level GLMs
run_higher_level    = False    # group-level analyses and statistics

# -----------------------
# Participants
# -----------------------
# Notes
# sub-211_ses-02 has no events file for RSA run 2 (need to adjust all first level functions)
# sub-107_ses-01 no T1 anatomical, have to overwrite in preprocessing FSF by hand to use session 2's T1
participants    = pd.read_csv(os.path.join(analysis_dir,'participants_debug.csv'), dtype=str) # open in textmate, not excel!
mri_subjects    = participants['mri_subjects']
subjects_group  = participants['subjects']
sessions        = ['ses-01','ses-02']
# the following is for indexing any missing files
t1              = ['t1_1','t1_2']
fieldmap        = ['fieldmap1','fieldmap2']
loc_letters     = ['loc-letters1','loc-letters2']
loc_colors      = ['loc-colors1','loc-colors2']
rsa             = [ ['rsa1_1','rsa2_1','rsa3_1','rsa4_1'],
                    ['rsa1_2','rsa2_2','rsa3_2','rsa4_2'] ] 

# -----------------------
# Run preprocessing class
# -----------------------
# Values for unwarping in FEAT, based on BOLD data
EPI_TE      = 0.0396*1000       # echo time EPI in ms
EPI_EECHO   = 0.000580009*1000  # effective echo spacing EPI in ms
TDIFF_ECHO  = 0.00246*1000      # difference between echo times of field map magnitude images [0.0047,0.00716] in ms
TR          = 1.5       # repetition time in seconds
FWHM        = 3         # smoothing kernel in mm (localizers only)

if run_preprocessing:
    for s,subject in enumerate(mri_subjects): # go in order of mri
        ######################################################################
        # check if field map exists in both sessions (identical preprocessing for each session)
        if int(participants[fieldmap[0]][s]) and int(participants[fieldmap[1]][s]):
            UNWARP = 1
        else: # otherwise, turn off unwarping for both sessions current subject
            UNWARP = 0
        ######################################################################
        for ss,session in enumerate(sessions): # loop sessions
            ###################################################################
            # check if T1 exists for current session. If not, use other session
            if int(participants[t1[ss]][s]):
                t1_sess = session
            else:
                t1_sess = sessions[~ss]
            
            ######################################################################
            # initialize class
            preprocess = letter_color_preprocessing.preprocess_class(
                subject         = subjects_group[s], # experiment subject number
                mri_subject     = subject.zfill(3),  # mri subject number, leading zeros
                session         = session,
                analysis_dir    = analysis_dir,
                source_dir      = source_dir,
                raw_dir         = raw_dir,
                deriv_dir       = deriv_dir,
                mask_dir        = mask_dir,
                template_dir    = template_dir,
                EPI_TE          = EPI_TE,
                EPI_EECHO       = EPI_EECHO,
                TDIFF_ECHO      = TDIFF_ECHO,
                UNWARP          = UNWARP,
                FWHM            = FWHM,
                t1_sess         = t1_sess
                )            
            # preprocess.dicom2bids()               # convert DICOMS from scanner to nifti in bids format
            # preprocess.housekeeping()             # copies events and physio files (rename event file names to be bids compliant and same b/t mri & behavior)
            # preprocess.raw_copy()                 # copy from bids_raw directory into derivaties to prevent overwriting
            # preprocess.sbref_raw_copy()           # hack, can DELETE later once finished
            # preprocess.bet_brains_T1()            # brain extraction T1 always check visually!
            # preprocess.bet_brains_fmap()          # B0 unwarping needs 'tight' brain extracted magnitude images of field map, better too small than too big!
            # preprocess.prepare_fmap()             # prepares the field map image in radians/sec
            
            #### Everything below here writes commands to a batch job ####
            
            # preprocess.preprocess_fsf()                   # generate FSF file for preprocessing in FEAT (run from command line - batch)
            
            #### AFTER FEAT FINISHES run:
            # preprocess.register_native2native_target()    # transform all runs to native_target, compute matrices but don't apply yet to bold time series
            # preprocess.invert_registrations()             # inverse and concatenate all necessary registration transformations
            # preprocess.register_ventricle2native()        # register the 4th ventricle (MNI space) to native space
            
            #### AFTER REGISTRATION run:
            # preprocess.ventricle_regressor()                # create a 4th ventricle regressor (NIFTI)
            # preprocess.motion_regressors()                  # create motion regressors (NIFTI)
            # preprocess.physiological_noise_regressors()     # creates the RETROICOR regressors from physiological recordings
            # preprocess.nuisance_regressor_list()            # create a list of all the nuisance regressors for the 1st level analysis
            
            preprocess.transform_native2native_target()    # transform all bold time series runs to native_target
            
            shell()
            
            
# -----------------------
# Run first-level class
# -----------------------   
if run_first_level:
    for s,subject in enumerate(subjects_group):
        # sessions are looped depending on task
        first_level = letter_color_first_level.first_level_class(
            subject         = subject,
            analysis_dir    = analysis_dir,
            deriv_dir       = deriv_dir,
            mask_dir        = mask_dir,
            template_dir    = template_dir,
            timing_files_dir = timing_files_dir,
            TR              = TR, # repitition time in seconds
            )
        first_level.loc_combine_epi('colors')             # concantenate both runs of localizer to perform 1 GLM
        # first_level.loc_combine_timing_files('colors')    # timing files for color localizer GLM
        # first_level.loc_nuisance_regressors('colors')     # concatenate motion parameters from preprocessing, also outputs cols of 1s for each blocks' mean
        # first_level.loc_fsf('colors')                     # generates the first level FSF for the localizers
        # 
        # first_level.loc_combine_epi('letters')            # concantenate both runs of localizer to perform 1 GLM
        # first_level.loc_combine_timing_files('letters')   # timing files for color localizer GLM
        # first_level.loc_nuisance_regressors('letters')    # concatenate motion parameters from preprocessing, also outputs cols of 1s for each blocks' mean
        # first_level.loc_fsf('letters')                    # generates the first level FSF for the localizers
        #
        # first_level.rsa_combine_epi()                     # concatenate EPI data for the 4 runs of the RSA task
        # first_level.rsa_combine_events()                  # concatenate events files for the 4 runs of the RSA task
        # first_level.rsa_nuisance_regressors()             # motion parameters, run means, and oddball trials as nuisance
        #
        # first_level.rsa_timing_files_letters()            # each letter in it's color and black: trained/untrained vs. color/black
        # first_level.rsa_letters_fsf()                     # generates the first level FSF for the RSA design
        
        # first_level.rsa_timing_files_2x2()                # simple 2x2 design: trained/untrained vs. color/black
        # first_level.rsa_2x2_fsf()                         # generates the first level FSF for the 2x2 design
                
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
