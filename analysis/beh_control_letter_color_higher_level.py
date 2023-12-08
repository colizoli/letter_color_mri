#!/usr/bin/env python
# encoding: utf-8

"""
higher_level.py
Created by O.Colizoli, June-2019
Last update: 11-06-2019
Python version 2.7

The following packages need to be installed because they are called from the command line, but not imported:
fsl

"""

import os, subprocess, sys
import shutil as sh
import pandas as pd
import numpy as np
import scipy as sp
from scipy import stats
import itertools 
from colorspace.colorlib import hexcols

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from IPython import embed as shell # for Oly's debugging only

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set(style='ticks', font='Arial', font_scale=1, rc={
    'axes.linewidth': 1,
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'text.color': 'Black',
    'axes.labelcolor':'Black',
    'xtick.color':'Black',
    'ytick.color':'Black',} )
sns.plotting_context()

############################################
# PLOT SIZES: (cols,rows)
# a single plot, 1 row, 1 col (2,2)
# 1 row, 2 cols (2*2,2*1)
# 2 rows, 2 cols (2*2,2*2)
# 2 rows, 3 cols (2*3,2*2)
# 1 row, 4 cols (2*4,2*1)
# Nsubjects rows, 2 cols (2*2,Nsubjects*2)


class higher_level_class(object):
    def __init__(self, subjects, sessions, analysis_dir, deriv_dir, participants):        
        self.subjects       = subjects
        self.sessions       = sessions
        self.analysis_dir   = str(analysis_dir)
        self.deriv_dir      = str(deriv_dir)
        self.participants   = participants # dataframe participants file
            
        self.higher_level_dir = os.path.join(self.deriv_dir,'higher_level')
        if not os.path.isdir(self.higher_level_dir):
            os.mkdir(self.higher_level_dir)

        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-questionnaires')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-questionnaires'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-iconic_memory')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-iconic_memory'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-consistency')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-consistency'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-choose_pairs')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-choose_pairs'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-stroop')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-stroop'))
            
        self.figure_dir =  os.path.join(self.higher_level_dir,'figures')
        if not os.path.isdir(self.figure_dir):
            os.mkdir(self.figure_dir)
        
        self.letters = [
            'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
        ]
        
        # Letter sets, from experiment/general_parameters.py
        self.letter_sets_lower = [
            ['b','c','d','e','g','n','o','p','q','r','w','x','y'], # 0 False 
            ['a','f','h','i','j','k','l','m','s','t','u','v','z']  # 1 True 
        ]
        
    
    def housekeeping(self,):
        # change folders and filenames from ses-1 to ses-01 
        
        bids_session = ['ses-01', 'ses-02']
                
        ###################
        # RENAME folder' 'ses-1' to 'ses-01' 
        ###################
        for s,subject in enumerate(self.subjects):
            try:
                for ss,session in enumerate(['ses-1', 'ses-2']):
                    dir_path = os.path.join(self.deriv_dir, 'sub-{}'.format(subject), session)
                    dir_new_path = os.path.join(self.deriv_dir, 'sub-{}'.format(subject), bids_session[ss])
                    os.rename(dir_path, os.path.join(dir_path, dir_new_path)) # old,new
                    print('old={} , new={}'.format(f,f_new))
            except:
                pass
        
        ###################
        # RENAME 'behav' 'ses-1' to 'ses-01' in behav logfile names
        ###################
        for s,subject in enumerate(self.subjects):
            try:
                for ss,session in enumerate(['ses-1', 'ses-2']):
                    dir_path = os.path.join(self.deriv_dir, 'sub-{}'.format(subject), bids_session[ss], 'behav')
                    for f in os.listdir(dir_path):                    
                        f_new = f.replace(session, bids_session[ss])
                        os.rename(os.path.join(dir_path, f),os.path.join(dir_path, f_new)) # old,new
                        print('old={} , new={}'.format(f,f_new))
            except:
                pass
        
        
    def dataframe_qualia(self,):    
        # calculate the PA score from the Reading Questionnaire
        # see PAscoring 2020.xls for details
        # PA score = projector - associator
        
        flip_scores = ['7', '15', '5', '6', '8', '16'] # flip the response direction of the Likert scale from 5 to 1 for these questions
        associator = ['7', '15', '5', '6', '8', '16', '10', '3', '14', '12', '18', '1'] # include the "not" questions
        projector = ['7', '15', '5', '6', '8', '16', '17', '2', '9', '11', '4', '13']
        
        ## use original scoring:
        # associator = [ '10', '3', '14', '12', '18', '1']
        # projector = ['17', '2', '9', '11', '4', '13']

        DF = pd.read_csv(os.path.join(self.higher_level_dir, 'task-questionnaires', 'participants_qualia.tsv'))
        DF = DF.loc[:, ~DF.columns.str.contains('^Unnamed')] # remove all unnamed columns
        
        DF2 = DF.copy()
        DF2[flip_scores] = np.max(DF[flip_scores]) + 1 - DF[flip_scores] # flip the likert scale
        
        DF2['associator'] = np.mean(DF[associator],axis=1)
        DF2['projector'] = np.mean(DF[projector],axis=1)
        DF['PA_score'] = DF2['projector'] - DF2['associator']
        
        ########## ADD COLUMN FOR GROUP ##########
        group = [
            # updating
            (DF['subject'] < 400), # group 1 300s
            (DF['subject'] >= 400) # group 2 400s
            ]
        values = [1,2]
        DF['group'] = np.select(group, values)
        
        DF.to_csv(os.path.join(self.higher_level_dir, 'task-questionnaires','participants_qualia.tsv'), sep='\t')
        print('success: dataframe_qualia')
    
    
    def plot_qualia(self,):
        # plots histograms of the PA scores, whole sample and per group
        # plots means of the PA scores, whole sample and per group
        
        ######### HISTOGRAM #########
        DF = pd.read_csv(os.path.join(self.higher_level_dir, 'task-questionnaires', 'participants_qualia.tsv'), sep='\t')
        DF.dropna(inplace=True)
        fig = plt.figure(figsize=(6,4))
        
        ######### whole sample #########
        ax = fig.add_subplot(231) 
        PA = np.array(DF['PA_score'])
        ax.hist(PA,bins=10,rwidth=2)
        ax.set_xlabel('PA score')
        ax.set_ylabel('Frequency')
        ax.set_title('N={}'.format(len(PA)))
        # ax.axvline(0, lw=1, alpha=0.3, color = 'k') # Add vertical line at t=0
        ax.set_xlim([-1.5,1.5]) # centered around 0
        ax.set_ylim([0,20]) # centered around 0
        
        ######### group 1 #########
        ax = fig.add_subplot(232) 
        PA1 = np.array(DF[DF['group'] == 1]['PA_score'])
        
        n,bins,patches = ax.hist(PA1,bins=10,rwidth=2)
        
        ax.set_xlabel('PA score')
        ax.set_ylabel('Frequency')
        ax.set_title('Group 1 N={}'.format(len(PA1)))
        # ax.axvline(0, lw=1, alpha=0.3, color = 'k') # Add vertical line at t=0
        ax.set_xlim([-1.5,1.5]) # centered around 0
        ax.set_ylim([0,20]) # centered around 0
        
        ######### group 2 #########
        ax = fig.add_subplot(233) 
        PA2 = np.array(DF[DF['group'] == 2]['PA_score'])
        # make sure to use the same bins as group 1
        ax.hist(PA2,bins=bins,rwidth=2)
        
        ax.set_xlabel('PA score')
        ax.set_ylabel('Frequency')
        ax.set_title('Group 2 N={}'.format(len(PA2)))
        # ax.axvline(0, lw=1, alpha=0.3, color = 'k') # Add vertical line at t=0
        ax.set_xlim([-1.5,1.5]) # centered around 0
        ax.set_ylim([0,20]) # centered around 0
        
        ######### MEANS #########
        ax = fig.add_subplot(234) 

        xlabels = ['All','Group 1','Group2']
        means = [np.mean(PA),np.mean(PA1),np.mean(PA2)]
        sems = [stats.sem(PA),stats.sem(PA1),stats.sem(PA2)]
        xind = np.arange(len(means))

        ax.bar(xind, means, yerr=sems, color='green',alpha=0.5)
        
        ax.set_xticklabels(xlabels)
        ax.set_ylabel('PA score')
        # stats.ttest_ind(PA1,PA2)
        
        ## whole figure
        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        fig.savefig(os.path.join(self.figure_dir,'histogram_PA_scores.pdf'))
        print('success: plot_qualia')
        
        
    def dataframe_choose_pairs(self,):
        # put all of group 1 and group 2's color choices together in a single dataframe
        DFOUT = pd.DataFrame()
        for s,subject in enumerate(self.subjects):
            if int(subject) < 400: # group 1 get colors file
                this_df = pd.read_csv(os.path.join(self.deriv_dir,'sub-{}'.format(subject),'sub-{}_colors.tsv'.format(subject)),sep='\t')
            else: # group 2 get preference file
                files = [f for f in os.listdir(os.path.join(self.deriv_dir, 'sub-{}'.format(subject), 'ses-01', 'behav')) if 'prefs' in f]
                if len(files) == 1:
                    this_df_path = os.path.join(self.deriv_dir, 'sub-{}'.format(subject), 'ses-01', 'behav', files[0])
                else:
                    print("Error check Preference logfiles sub-{}".format(subject))
                this_df = pd.read_csv(this_df_path, sep='\t')
            this_df = this_df.loc[:, ~this_df.columns.str.contains('^Unnamed')] # remove all unnamed columns
            DFOUT = pd.concat([DFOUT,this_df],axis=0) #concant on rows\
        
        ########## ADD COLUMN FOR GROUP ##########
        group = [
            # updating
            (DFOUT['subject'] < 400), # group 1
            (DFOUT['subject'] >= 400) # group 2
            ]
        values = [1,2]
        DFOUT['group'] = np.select(group, values)
        
        ########## ADD UNIQUE IDENTIFIERS FOR EACH COLOR ##########
        # copied from first_level.rsa_combine_events()
        rgb_codes = [
            (DFOUT['r'] == 0) & (DFOUT['g'] == 163) & (DFOUT['b'] == 228), # light_blue
            (DFOUT['r'] == 161) & (DFOUT['g'] == 199) & (DFOUT['b'] == 70), # lime_green
            (DFOUT['r'] == 183) & (DFOUT['g'] == 61) & (DFOUT['b'] == 160), # magenta
            (DFOUT['r'] == 181) & (DFOUT['g'] == 44) & (DFOUT['b'] == 67), # dark_red
            (DFOUT['r'] == 16) & (DFOUT['g'] == 114) & (DFOUT['b'] == 86), # dark_green
            (DFOUT['r'] == 237) & (DFOUT['g'] == 114) & (DFOUT['b'] == 162), # pink
            (DFOUT['r'] == 58) & (DFOUT['g'] == 175) & (DFOUT['b'] == 75), # green
            (DFOUT['r'] == 248) & (DFOUT['g'] == 154) & (DFOUT['b'] == 28), # light_orange
            (DFOUT['r'] == 109) & (DFOUT['g'] == 57) & (DFOUT['b'] == 142), # purple
            (DFOUT['r'] == 239) & (DFOUT['g'] == 79) & (DFOUT['b'] == 41), # orange
            (DFOUT['r'] == 49) & (DFOUT['g'] == 60) & (DFOUT['b'] == 163), # blue
            (DFOUT['r'] == 255) & (DFOUT['g'] == 211) & (DFOUT['b'] == 0), # yellow
            (DFOUT['r'] == 9) & (DFOUT['g'] == 181) & (DFOUT['b'] == 172) # teal
        ]
        
        hex_codes = [
            '#00A3E4',
            '#A1C746',
            '#B73DA0',
            '#B52C43',
            '#107256',
            '#ED72A2',
            '#3AAF4B',
            '#F89A1C',
            '#6D398E',
            '#EF4F29',
            '#313CA3',
            '#FFD300',
            '#09B5AC'
        ]
        
        color_names = [
            'light_blue',
            'lime_green',
            'magenta',
            'dark_red',
            'dark_green',
            'pink',
            'green',
            'light_orange',
            'purple',
            'orange',
            'blue',
            'yellow',
            'teal'
        ]
        DFOUT['hex_codes'] = np.select(rgb_codes, hex_codes)
        DFOUT['color_name'] = np.select(rgb_codes, color_names)
        
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-choose_pairs','task-choose_pairs_subjects.tsv'),sep='\t')
        print('success: dataframe_choose_pairs')
        
        
    def plot_choose_pairs(self,):
        # pie charts for frequency of letter-color pairs
        # group by letter and group1 vs. group2 (preference)
        DF = pd.read_csv(os.path.join(self.higher_level_dir, 'task-choose_pairs', 'task-choose_pairs_subjects.tsv'),sep='\t')
        DF = DF.loc[:, ~DF.columns.str.contains('^Unnamed')] # remove all unnamed columns
        
        G = pd.DataFrame(DF.groupby(['letter','group','hex_codes','color_name'])['subject'].count())
        G.reset_index(inplace=True)

        fig = plt.figure(figsize=(24,12))
        # loop through letters for each group
        counter = 1 # subplots
        ########## group 1 ##########
        for letter_set in self.letter_sets_lower:
            G1 = G[G['group']==1]
    
            for L in letter_set:
                ax = fig.add_subplot(5,13,counter)
                G1L = G1[G1['letter']==L]
                ax.pie(np.array(G1L['subject']),startangle=90, colors=np.array(G1L['hex_codes']))
                ax.set_title('Group 1, {}'.format(L))
                counter += 1
        
        ########## group 2 ##########
            G2 = G[G['group']==2]
            for L in letter_set:
                ax = fig.add_subplot(5,13,counter) # shift over and down
                G2L = G2[G2['letter']==L]
                ax.pie(np.array(G2L['subject']),startangle=90, colors=np.array(G2L['hex_codes']))
                ax.set_title('Group 2, {}'.format(L))
                counter += 1
        # last subplot show all colors
        ax = fig.add_subplot(5,13,counter) 
        all_colors = np.unique(np.array(G2L['hex_codes']))
        ax.pie(np.repeat(1,len(all_colors)),startangle=90, colors=all_colors)
        ax.set_title('color choices')
        
        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        fig.savefig(os.path.join(self.figure_dir,'frequency_choose_pairs.pdf'))
        
        print('success: plot_choose_pairs')
        
        
    def dataframe_subjects_iconic_memory(self,):
        # concantenates all behavioral logfiles of the iconic memory task
        # saves as large TSV file in higher_level/iconic_memory directory
        # counts missing, then drops missings trials
        # adds column for group membership
        
        DFOUT = pd.DataFrame()
        for s, subject in enumerate(self.subjects):
            for ss, session in enumerate(self.sessions):
                files = [f for f in os.listdir(os.path.join(self.deriv_dir, 'sub-{}'.format(subject), session, 'behav')) if 'task-iconic_behav' in f]
                if len(files) == 1:
                    this_df_path = os.path.join(self.deriv_dir, 'sub-{}'.format(subject), session, 'behav', files[0])
                else:
                    print("Error check Iconic Memory logfiles:  sub-{}".format(subject))
                this_df = pd.read_csv(this_df_path, dtype={'subject':np.int32, 'session':np.int32, 'letter_set':np.int32, 'trained':np.int32, 'trial_number':np.int32,
                'block_trial_number':np.int32,'target_position':np.int32,'correct':np.int32},sep='\t')
                
                DFOUT = pd.concat([DFOUT,this_df],axis=0) #concant on rows
        DFOUT = DFOUT.loc[:, ~DFOUT.columns.str.contains('^Unnamed')] # remove all unnamed columns
        
        ########## COUNT MISSING TRIALS PER SUBJECT ##########
        MISSING = DFOUT.groupby(['subject','session']).count().rsub(DFOUT.groupby(['subject','session']).size(), axis=0)['RT']
        MISSING = pd.DataFrame(MISSING)
        MISSING.to_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_missing.tsv'),sep='\t')
        
        ########## DROP MISSING TRIALS ##########
        DFOUT.dropna(subset=['RT'],inplace=True)
    
        ########## ADD COLUMN FOR GROUP ##########
        group = [
            # updating
            (DFOUT['subject'] < 400), # group 1
            (DFOUT['subject'] >= 400) # group 2
            ]
        values = [1,2]
        DFOUT['group'] = np.select(group, values)

        # resave log file with new columns in higher_level folder
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_subjects.tsv'),sep='\t')
        print('success: dataframe_subjects_iconic_memory')
      
      
    def dataframe_anova_iconic_memory(self,):
        # create dataframe for JASP input: mixed ANOVA
        # group (1 vs. 2), session (1 vs. 2), letter condition (trained vs. untrained)
        # columns: ses-01_untrained	ses-01_trained	ses-02_untrained ses-02_trained
        
        DFIN = pd.read_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_subjects.tsv'),sep='\t')
        
        # groupby
        G = DFIN.groupby(['group','session','trained','subject'])['correct'].mean()
        # unstack: subjects as rows, cols = session x letter condition and group membership
        G0 = G.unstack(level=1)
        G1 = G0.unstack(level=1)
        
        # drop first two rows and replace cols with: ses-01_untrained	ses-01_trained	ses-02_untrained ses-02_trained	
        G1.columns = G1.columns.droplevel()
        G1.columns = ['session1_untrained', 'session1_trained',	'session2_untrained', 'session2_trained']
        G1.to_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_anova.tsv'),sep='\t')
        print('success: dataframe_anova_iconic_memory')
    
    
    def plot_anova_iconic_memory(self,):
        # plot the bar graphs for the 2x2 ANOVA (collapse over groups)
        
        DF = pd.read_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_anova.tsv'),sep='\t')
        DF = DF.multiply(100) # for percentage
        
        labels = ['Pre', 'Post']
        untrained_means = [DF.mean()[2],DF.mean()[4]]
        trained_means = [DF.mean()[3],DF.mean()[5]]
            
        untrained_sems = [stats.sem(DF)[2], stats.sem(DF)[4]]
        trained_sems = [stats.sem(DF)[3], stats.sem(DF)[5]]

        x = np.arange(len(labels))  # the label locations
        width = 0.45  # the width of the bars

        fig = plt.figure(figsize=(2,2))
        ax = fig.add_subplot(111)
                
        rects1 = ax.bar(x - width/2, untrained_means, yerr=untrained_sems, width=width, label='Untrained',color='black',alpha=0.3,edgecolor=None,linewidth=0)
        rects2 = ax.bar(x + width/2, trained_means, yerr=trained_sems, width=width, label='Trained',color='black',alpha=0.6,edgecolor=None,linewidth=0)
        # chance level
        # ax.axhline(0.125, lw=1, alpha=0.3, color = 'k') # Add horizontal line at t=0
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylim([30,50])
        ax.set_ylabel('Accuracy')
        ax.set_title('Iconic Memory')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(loc='upper left')

        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        fig.savefig(os.path.join(self.figure_dir,'anova_iconic_memory.pdf'))
        print('success: plot_anova_iconic_memory')
        
        
    def dataframe_subjects_cieluv(self,):
        # concantenates all behavioral logfiles of the consistency task
        # converts RGB to CIELUV space (save as separate data frame in higher_level/consistency directory)
        
        ################################################################################################
        # converts RGB to CIELUV space (save as separate data frame in higher_level/consistency directory)
        DFOUT = pd.DataFrame() 
        for s,subject in enumerate(self.subjects):
            # check if missing file either session
            for ss,session in enumerate(self.sessions):
                print(subject)
                print(session)
                
                files = [f for f in os.listdir(os.path.join(self.deriv_dir, 'sub-{}'.format(subject), session, 'behav')) if 'task-consistency' in f]
                if len(files) == 1:
                    this_df_path = os.path.join(self.deriv_dir, 'sub-{}'.format(subject), session, 'behav', files[0])
                else:
                    print("Error check Consistency logfiles:  sub-{}".format(subject))
                    
                this_df_path = os.path.join(this_df_path)

                try: # have for comma vs. tab separated
                    this_df = pd.read_csv(this_df_path,header=0)
                except:
                    this_df = pd.read_csv(this_df_path,header=0,sep='\t')
                if len(this_df.columns.values) < 2:
                    this_df = pd.read_csv(this_df_path,header=0,sep='\t')
                    
                ########## DROP MISSING LETTERS ##########
                this_df.dropna(subset=['hex'],inplace=True)
                ########## COLOR CONVERSION ##########
                c = hexcols(this_df['hex'])
                # c.set_whitepoint('#FFFFFF')
                c.to("CIELUV") # covert to CIELUV space, happens in place
                LUV = c.get() # dictionary object for each dimension
                L = list(LUV.values())[0] # CIELUV values
                U = list(LUV.values())[1]
                V = list(LUV.values())[2]
                this_df[list(LUV)[0]] = L # store as separate columns
                this_df[list(LUV)[1]] = U
                this_df[list(LUV)[2]] = V
                this_df['session'] = np.repeat(ss+1,len(this_df))
                this_df = this_df.loc[:, ~this_df.columns.str.contains('^Unnamed')] # remove all unnamed columns
                
                ########## ADD COLUMN FOR GROUP ##########
                group = [
                    # updating
                    (this_df['subject'] < 400), # group 1
                    (this_df['subject'] >= 400) # group 2
                    ]
                values = [1,2]
                this_df['group'] = np.select(group, values)
                
                DFOUT = pd.concat([DFOUT,this_df[['subject','session','group','choice','letter','hex','rgb','L','U','V']]],axis=0) #concant on rows
            DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_cieluv_subjects.tsv'),sep='\t') # save
        print('success: dataframe_subjects_cieluv')
        
        
    def dataframe_subjects_consistency(self,):        
        # calculates Euclidean distance between choice 1 and 2 for each letter, session, subject. 
        # Save dataframe higher_level/consistency directory)
        # adds column for group membership

        CHOICES = pd.read_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_cieluv_subjects.tsv'),sep='\t') # save
        DFOUT = pd.DataFrame(columns=['subject','session','group','letter','trained','distance']) 
        counter = 0
        for s,subject in enumerate(np.unique(CHOICES['subject'])):         
            
            # load colors file to match letters to training conditions
            this_colors = pd.read_csv(os.path.join(self.deriv_dir,'sub-{}'.format(subject),'sub-{}_colors.tsv'.format(subject)),sep='\t')
               
            for ss,session in enumerate(np.unique(CHOICES['session'])):
                # select current data
                    this_df1 = CHOICES[(CHOICES['subject']==subject) & (CHOICES['session']==session) & (CHOICES['choice']==1)] # select choice 1
                    this_df2 = CHOICES[(CHOICES['subject']==subject) & (CHOICES['session']==session) & (CHOICES['choice']==2)] # select choice 2

                    # loop through letters, check whether all letters are included
                    for L in list(set(np.array(this_df1['letter'])).intersection(np.array(this_df2['letter']))):
                        
                        ########## DISTANCE ##########
                        this_letter1 = this_df1[(this_df1['letter']==L)]
                        this_letter2 = this_df2[(this_df2['letter']==L)]
                        a = np.array(this_letter1[['L','U','V']])
                        b = np.array(this_letter2[['L','U','V']])
                        EDIST = np.sqrt(np.sum(np.square(a-b)))
                        
                        ########## TRAINED LETTER? ##########
                        trained = len(list(set(L).intersection(np.array(this_colors['letter'])))) > 0
                        
                        ########## ADD ROW TO OUTPUT DATAFRAME ##########
                        DFOUT.loc[counter] = [
                            subject,            # subject number
                            session,            # session
                            np.unique(this_df1['group'])[0],   # group
                            L,                  # letter 
                            trained,            # trained letter?
                            EDIST,              # Euclidean distance in CIELUV space
                        ]          
                        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_cieluv_letters.tsv'),sep='\t')
                        counter += 1        
        print('success: dataframe_subjects_consistency')
    
    
    def dataframe_anova_consistency(self,):
        # create dataframe for JASP input: mixed ANOVA
        # group (1 vs. 2), session (1 vs. 2), letter condition (trained vs. untrained)
        
        DFIN = pd.read_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_cieluv_letters.tsv'),sep='\t')
        
        # groupby
        G = DFIN.groupby(['group','session','trained','subject'])['distance'].mean()
        # unstack: subjects as rows, cols = session x letter condition and group membership
        G0 = G.unstack(level=1)
        G1 = G0.unstack(level=1)
        
        # drop first two rows and replace cols with: ses-01_untrained	ses-01_trained	ses-02_untrained ses-02_trained	
        G1.columns = G1.columns.droplevel()
        G1.columns = ['session1_untrained', 'session1_trained',	'session2_untrained', 'session2_trained']
        
        G1.to_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_anova.tsv'),sep='\t')
        print('success: dataframe_anova_consistency')
    
    
    def plot_anova_consistency(self,):
        # plot the bar graphs for the 2x2 ANOVA (collapse over groups)
        
        DF = pd.read_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_anova.tsv'),sep='\t')

        DF.dropna(inplace=True)
        
        labels = ['Pre', 'Post']
        untrained_means = [DF.mean()[2],DF.mean()[4]]
        trained_means = [DF.mean()[3],DF.mean()[5]]
            
        untrained_sems = [stats.sem(DF)[2],stats.sem(DF)[4]]
        trained_sems = [stats.sem(DF)[3],stats.sem(DF)[5]]

        x = np.arange(len(labels))  # the label locations
        width = 0.45  # the width of the bars

        fig = plt.figure(figsize=(2,2))
        ax = fig.add_subplot(111)
        
        rects1 = ax.bar(x - width/2, untrained_means, yerr=untrained_sems, width=width, label='Untrained',color='black',alpha=0.3,edgecolor=None,linewidth=0)
        rects2 = ax.bar(x + width/2, trained_means, yerr=trained_sems, width=width, label='Trained',color='black',alpha=0.6,edgecolor=None,linewidth=0)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        # ax.set_ylim([0.3,0.5])
        ax.set_ylabel('Distance')
        ax.set_title('Consistency')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(loc='lower left')

        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        fig.savefig(os.path.join(self.figure_dir,'anova_consistency.pdf'))
        print('success: plot_anova_consistency')
        
        
    def correlation_consistency_iconic_memory(self,):
        # test the correlation in the letter conditions between consistency scores and iconic memory performance 
        
        IM = pd.read_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_anova.tsv'),sep='\t')
        CON = pd.read_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_anova.tsv'),sep='\t')
        
        # remove missing subjects        
        # missing = list(set(np.array(IM['subject'])).difference(np.array(CON['subject'])))
        # IM = IM[np.array(IM['subject'])!=missing]
        IM = IM.multiply(100) # for percentage 
        # CON = CON[np.array(CON['subject'])!=missing]
        #
        fig = plt.figure(figsize=(8,4))
        
        ########################
        # subplot 1 trained
        ax = fig.add_subplot(121)
        # PRE TRAINING
        x = np.array(CON['session1_trained'])
        y = np.array(IM['session1_trained'])
        r,pval = stats.pearsonr(x,y)
        ax.plot(x, y, 'o',color='purple',alpha=0.5,label='pre r={},p={}'.format(np.round(r,2),np.round(pval,3)))
        m, b = np.polyfit(x, y, 1)
        ax.plot(x, m*x+b,color='purple',alpha=0.5)
        # POST TRAINING
        x = np.array(CON['session2_trained'])
        y = np.array(IM['session2_trained'])
        r,pval = stats.pearsonr(x,y)
        ax.plot(x, y, 'o',color='purple',alpha=1, label='post r={},p={}'.format(np.round(r,2),np.round(pval,3)))
        m, b = np.polyfit(x, y, 1)
        ax.plot(x, m*x+b,color='purple',alpha=1)
        
        ax.set_ylim([0,75])
        ax.set_xlim([10,160])
        ax.set_xlabel('Consistency color distance')
        ax.set_ylabel('Iconic memory accuracy')
        ax.set_title('Trained letters N={}'.format(len(x)))
        ax.legend(loc='lower left')
                
        ########################
        # subplot 2 untrained
        ax = fig.add_subplot(122)

        x = np.array(CON['session1_untrained'])
        y = np.array(IM['session1_untrained'])
        r,pval = stats.pearsonr(x,y)
        ax.plot(x, y, 'o',color='purple',alpha=0.5,label='pre r={},p={}'.format(np.round(r,2),np.round(pval,3)))
        m, b = np.polyfit(x, y, 1)
        ax.plot(x, m*x+b,color='purple',alpha=0.5)
        
        x = np.array(CON['session2_untrained'])
        y = np.array(IM['session2_untrained'])
        r,pval = stats.pearsonr(x,y)
        ax.plot(x, y, 'o',color='purple',alpha=1,label='post r={},p={}'.format(np.round(r,2),np.round(pval,3)))
        m, b = np.polyfit(x, y, 1)
        ax.plot(x, m*x+b,color='purple',alpha=1)
        
        ax.set_ylim([0,75])
        ax.set_xlim([10,160])
        ax.set_xlabel('Consistency color distance')
        ax.set_ylabel('Iconic memory accuracy')
        ax.set_title('Untrained letters N={}'.format(len(x)))
        ax.legend(loc='lower left')
                
        # whole figure format
        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        fig.savefig(os.path.join(self.figure_dir,'correlation_consistency_iconic_memory.pdf'))
        print('success: correlation_consistency_iconic_memory')
        

    def dataframe_subjects_stroop(self,):
        # concantenates all behavioral logfiles of the stroop task
        # saves as large TSV file in higher_level/stroop directory
        # counts missing, then drops missings trials
        # adds column for group membership
        
        DFOUT = pd.DataFrame()
        for s, subject in enumerate(self.subjects):
            for ss, session in enumerate(self.sessions):
                files = [f for f in os.listdir(os.path.join(self.deriv_dir, 'sub-{}'.format(subject), session, 'behav')) if 'stroop' in f]
                if len(files) == 1:
                    this_df_path = os.path.join(self.deriv_dir, 'sub-{}'.format(subject), session, 'behav', files[0])
                else:
                    print("Error check Stroop logfiles:  sub-{}".format(subject))
                print(subject)
                print(session)
                try:
                    this_df = pd.read_csv(this_df_path, sep='\t')
                except:
                    continue
            
                DFOUT = pd.concat([DFOUT,this_df],axis=0) #concant on rows

        DFOUT = DFOUT.loc[:, ~DFOUT.columns.str.contains('^Unnamed')] # remove all unnamed columns
        
        ########## COUNT MISSING TRIALS PER SUBJECT ##########
        MISSING = DFOUT.groupby(['subject','session']).count().rsub(DFOUT.groupby(['subject','session']).size(), axis=0)['RT']
        MISSING = pd.DataFrame(MISSING)
        MISSING.to_csv(os.path.join(self.higher_level_dir,'task-stroop','task-stroop_missing.tsv'),sep='\t')
        
        ########## DROP MISSING TRIALS ##########
        DFOUT.dropna(subset=['RT'],inplace=True)
    
        ########## ADD COLUMN FOR GROUP ##########
        group = [
            # updating
            (DFOUT['subject'] < 400), # group 1
            (DFOUT['subject'] >= 400) # group 2
            ]
        values = [1,2]
        DFOUT['group'] = np.select(group, values)

        # resave log file with new columns in higher_level folder
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-stroop','task-stroop_subjects.tsv'),sep='\t')
        print('success: dataframe_subjects_stroop')


    def dataframe_anova_stroop(self,):
        # create dataframe for JASP input: mixed ANOVA
        # group (1 vs. 2), session (1 vs. 2), letter condition (trained vs. untrained)
        # columns: 'session1_incongruent', 'session1_congruent',	'session2_incongruent', 'session2_congruent'
        
        DFIN = pd.read_csv(os.path.join(self.higher_level_dir,'task-stroop','task-stroop_subjects.tsv'),sep='\t')
        
        # groupby
        G = DFIN.groupby(['group','session','congruent','subject'])['RT'].mean()
        # unstack: subjects as rows, cols = session x letter condition and group membership
        G0 = G.unstack(level=1)
        G1 = G0.unstack(level=1)
        
        # drop first two rows and replace cols with: ses-01_untrained	ses-01_trained	ses-02_untrained ses-02_trained	
        G1.columns = G1.columns.droplevel()
        G1.columns = ['session1_incongruent', 'session1_congruent',	'session2_incongruent', 'session2_congruent']
        G1.to_csv(os.path.join(self.higher_level_dir,'task-stroop','task-stroop_anova.tsv'),sep='\t')
        print('success: dataframe_anova_stroop')
    

    def plot_anova_stroop(self,):
        # plot the bar graphs for the 2x2 ANOVA (collapse over groups)
        
        DF = pd.read_csv(os.path.join(self.higher_level_dir,'task-stroop','task-stroop_anova.tsv'),sep='\t')
        
        DF.dropna(inplace=True)
        
        labels = ['Pre', 'Post']
        untrained_means = [DF.mean()[2],DF.mean()[4]]
        trained_means = [DF.mean()[3],DF.mean()[5]]
            
        untrained_sems = [stats.sem(DF)[2], stats.sem(DF)[4]]
        trained_sems = [stats.sem(DF)[3], stats.sem(DF)[5]]

        x = np.arange(len(labels))  # the label locations
        width = 0.45  # the width of the bars

        fig = plt.figure(figsize=(2,2))
        ax = fig.add_subplot(111)
                
        rects1 = ax.bar(x - width/2, untrained_means, yerr=untrained_sems, width=width, label='Incongruent',color='black',alpha=0.3,edgecolor=None,linewidth=0)
        rects2 = ax.bar(x + width/2, trained_means, yerr=trained_sems, width=width, label='Congruent',color='black',alpha=0.6,edgecolor=None,linewidth=0)
        # chance level
        # ax.axhline(0.125, lw=1, alpha=0.3, color = 'k') # Add horizontal line at t=0
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylim([0.75,0.9])
        ax.set_ylabel('RT')
        ax.set_title('Stroop')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(loc='upper left')

        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        fig.savefig(os.path.join(self.figure_dir,'anova_stroop.pdf'))
        print('success: plot_anova_stroop')
        