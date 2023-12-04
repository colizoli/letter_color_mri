#!/usr/bin/env python
# encoding: utf-8

"""
run_beh_control_letter_color_analysis.py
Created by O.Colizoli, 2023
Last update: 04-12-2023
Python version 3.6

Behavioral control study including a single-letter Stroop task

"""

# -----------------------
# Import functionality
# ----------------------- 
import os, subprocess, sys
import pandas as pd
import numpy as np
from IPython import embed as shell # for Oly's debugging only
# custom analysis scripts
import beh_control_letter_color_higher_level as letter_color_higher_level

# -----------------------
# Paths
# ----------------------- 
##########################
home_dir = os.path.dirname(os.getcwd()) # one level up from analysis folder
# home_dir = '/project/3018051.01/'
##########################
analysis_dir    = os.path.join(home_dir, 'analysis')        # scripts + configuration files for BIDS
source_dir      = os.path.join(home_dir, 'raw')             # as they were acquired
deriv_dir       = os.path.join(home_dir,'derivatives')      # Processed data

# -----------------------
# Levels (switch ON/OFF)
# ----------------------- 
run_higher_level    = True    # group-level analyses and statistics

# -----------------------
# Participants
# -----------------------
# Notes
# sub-211_ses-02 has no events file for RSA run 2 (need to adjust all first level functions)
# sub-107_ses-01 no T1 anatomical, have to overwrite in preprocessing FSF by hand to use session 2's T1
participants = pd.read_csv(os.path.join(analysis_dir,'participants_beh_control.csv'), dtype=str) # open in textmate, not excel!
subjects     = participants['subject']
sessions     = ['ses-01','ses-02']
                
# -----------------------
# Run higher-level class
# -----------------------
if run_higher_level:
    ### PARTICIPANT LIST FOR BEHAVIOR ONLY
    higher_level = letter_color_higher_level.higher_level_class(
            subjects     = subjects,
            sessions     = sessions,
            analysis_dir = analysis_dir,
            deriv_dir    = deriv_dir,
            participants = participants
            )
            
    # higher_level.housekeeping()                             # ses-1 to ses-01
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

    # higher_level.dataframe_subjects_stroop()
    # higher_level.dataframe_anova_stroop()
    higher_level.plot_anova_stroop()                 # plot ANOVA 2x2 (collapsed over groups)
    

