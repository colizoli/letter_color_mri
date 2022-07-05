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
import nibabel as nib
import pandas as pd
import numpy as np
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
    def __init__(self, subjects, sessions, analysis_dir, deriv_dir, mask_dir, template_dir, TR, participants):        
        self.subjects       = subjects
        self.sessions       = sessions
        self.analysis_dir   = str(analysis_dir)
        self.deriv_dir      = str(deriv_dir)
        self.mask_dir       = str(mask_dir)
        self.template_dir   = str(template_dir)
        self.TR             = str(TR)
        self.participants   = participants # dataframe participants file

        if not os.path.isdir(self.mask_dir):
            os.mkdir(self.mask_dir)
        if not os.path.isdir(self.template_dir):
            os.mkdir(self.template_dir)
            
        self.preprocess_dir = os.path.join(self.deriv_dir,'preprocessing')
        self.first_level_dir = os.path.join(self.deriv_dir,'first_level')
        
        self.higher_level_dir = os.path.join(self.deriv_dir,'higher_level')
        if not os.path.isdir(self.higher_level_dir):
            os.mkdir(self.higher_level_dir)
            
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-colors')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-colors'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-letters')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-letters'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-rsa')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-rsa'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-iconic_memory')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-iconic_memory'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-consistency')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-consistency'))
        if not os.path.isdir(os.path.join(self.higher_level_dir,'task-choose_pairs')):
            os.mkdir(os.path.join(self.higher_level_dir,'task-choose_pairs'))
            
        self.figure_dir =  os.path.join(self.higher_level_dir,'figures')
        if not os.path.isdir(self.figure_dir):
            os.mkdir(self.figure_dir)
            
        # brain atlases 
        # note: 'maxprob-thr0-2mm' files have labels in single volume, while the 'prob-2mm' files contain probabilities for each label as time axis
        self.mni        = os.path.join(self.mask_dir,'MNI152_T1_2mm.nii.gz')
        self.mni_brain  = os.path.join(self.mask_dir,'MNI152_T1_2mm_brain.nii.gz')
        self.mni_labels = os.path.join(self.mask_dir,'atlases','MNI','MNI-maxprob-thr0-2mm.nii.gz')
        self.ho_cortical_labels = os.path.join(self.mask_dir,'atlases','HarvardOxford','HarvardOxford-{}-maxprob-thr0-2mm.nii.gz'.format('cort'))
        self.ho_cortical_prob = os.path.join(self.mask_dir,'atlases','HarvardOxford','HarvardOxford-{}-prob-2mm.nii.gz'.format('cort'))
        
        # write unix commands to job to run in parallel
        self.higher_level_job_path = os.path.join(self.analysis_dir,'jobs','job_higher_level.txt')
        if not os.path.exists(self.higher_level_job_path):
            self.higher_level_job = open(self.higher_level_job_path, "w")
            self.higher_level_job.write("#!/bin/bash\n")
            self.higher_level_job.close()
        
        self.letters = [
            'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
        ]
        
        # Letter sets, from experiment/general_parameters.py
        self.letter_sets_lower = [
            ['b','c','d','e','g','n','o','p','q','r','w','x','y'], # 0 False 
            ['a','f','h','i','j','k','l','m','s','t','u','v','z']  # 1 True 
        ]
    
    def dataframe_choose_pairs(self,):
        # put all of group 1 and group 2's color choices together in a single dataframe
        DFOUT = pd.DataFrame()
        for s,subject in enumerate(self.subjects):
            if int(subject) < 200: # group 1 get colors file
                this_df = pd.read_csv(os.path.join(self.deriv_dir,'sub-{}'.format(subject),'sub-{}_colors.tsv'.format(subject)),sep='\t')
            else: # group 2 get preference file
                this_df = pd.read_csv(os.path.join(self.deriv_dir,'sub-{}'.format(subject),'ses-01','behav','sub-{}_prefs.tsv'.format(subject,'ses-01')),sep='\t')
            this_df = this_df.loc[:, ~this_df.columns.str.contains('^Unnamed')] # remove all unnamed columns
            DFOUT = pd.concat([DFOUT,this_df],axis=0) #concant on rows\
        
        ########## ADD COLUMN FOR GROUP ##########
        group = [
            # updating
            (DFOUT['subject'] < 200), # group 1
            (DFOUT['subject'] >= 200) # group 2
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
        DF = pd.read_csv(os.path.join(self.higher_level_dir,'task-choose_pairs','task-choose_pairs_subjects.tsv'),sep='\t')
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
        for s,subject in enumerate(self.subjects):
            for ss,session in enumerate(self.sessions):
                this_df_path = os.path.join(self.deriv_dir,'sub-{}'.format(subject),session,'behav','sub-{}_{}_task-iconic_behav.tsv'.format(subject,session))
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
            (DFOUT['subject'] < 200), # group 1
            (DFOUT['subject'] >= 200) # group 2
            ]
        values = [1,2]
        DFOUT['group'] = np.select(group, values)

        # resave log file with new columns in higher_level folder
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_subjects.tsv'),sep='\t')
        print('success: dataframe_subjects_iconic_memory')
      
    def dataframe_anova_iconic_memory(self,):
        # create dataframe for JASP input: mixed ANOVA
        # group (1 vs. 2), session (1 vs. 2), letter condition (trained vs. untrained)
        
        DFIN = pd.read_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_subjects.tsv'),sep='\t')
        
        # groupby
        G = DFIN.groupby(['group','session','trained','subject'])['correct'].mean()
        # unstack: subjects as rows, cols = session x letter condition and group membership
        G0 = G.unstack(level=1)
        G1 = G0.unstack(level=1)
        G1.to_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_anova.tsv'),sep='\t')
        print('success: dataframe_anova_iconic_memory')
    
    def plot_anova_iconic_memory(self,):
        # plot the bar graphs for the 2x2 ANOVA (collapse over groups)
        
        DF = pd.read_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_anova.tsv'),sep='\t')
        DF = DF.multiply(100) # for percentage
        
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
        print('success: dataframe_anova_iconic_memory')
        
    def dataframe_subjects_cieluv(self,):
        # concantenates all behavioral logfiles of the consistency task
        # converts RGB to CIELUV space (save as separate data frame in higher_level/consistency directory)
        
        ################################################################################################
        # converts RGB to CIELUV space (save as separate data frame in higher_level/consistency directory)
        DFOUT = pd.DataFrame() 
        for s,subject in enumerate(self.subjects):
            # check if missing file either session
            if int(self.participants['consistency1'][s]) and int(self.participants['consistency2'][s]): 
                for ss,session in enumerate(self.sessions):
                    print(subject)
                    print(session)
                    this_df_path = os.path.join(self.deriv_dir,'sub-{}'.format(subject),session,'behav','sub-{}_{}_task-consistency_events.tsv'.format(subject,session))

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
                        (this_df['subject'] < 200), # group 1
                        (this_df['subject'] >= 200) # group 2
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
        G1.to_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_anova.tsv'),sep='\t')
        print('success: dataframe_anova_consistency')
    
    def plot_anova_consistency(self,):
        # plot the bar graphs for the 2x2 ANOVA (collapse over groups)
        
        DF = pd.read_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_anova.tsv'),sep='\t')

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
        print('success: dataframe_anova_consistency')
        
    def correlation_consistency_iconic_memory(self,):
        # test the correlation in the letter conditions between consistency scores and iconic memory performance 
        
        IM = pd.read_csv(os.path.join(self.higher_level_dir,'task-iconic_memory','task-iconic_memory_anova.tsv'),sep='\t')
        CON = pd.read_csv(os.path.join(self.higher_level_dir,'task-consistency','task-consistency_anova.tsv'),sep='\t')
        
        # remove missing subjects        
        missing = list(set(np.array(IM['subject'])).difference(np.array(CON['subject'])))
        IM = IM[np.array(IM['subject'])!=missing]
        IM = IM.multiply(100) # for percentage 
        CON = CON[np.array(CON['subject'])!=missing]
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
        
    def localizers_randomise_input(self,task):
        # Localizers (letters, colors), cope1 was contrast of interest
        # For randomise input, we need the 3D cope1 images stacked for all subjects in the 4th dimension
        # i.e., MNI x #subjects
        
        if task == 'letters':
            mask = os.path.join(self.mask_dir,'occipital_fusiform_temporal_occipital_fusiform_cortex_LH.nii.gz') # search only within these voxels
        elif task == 'colors':
            mask =  os.path.join(self.mask_dir,'occipital_temporal_occipital_fusiform_cortex.nii.gz')
        
        # output path MNI x #subjects
        outFile = os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_cope1.nii.gz'.format(task))
        # Load MNI brain and save headers
        mni = nib.load(self.mni)
        # initialize empty 4D array
        COPE = np.zeros((mni.get_data().shape[0],mni.get_data().shape[1],mni.get_data().shape[2],len(self.subjects)))
        for s,subject in enumerate(self.subjects):
            C1_path = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}.feat'.format(task,subject),'stats','cope1.nii.gz')
            C1 = nib.load(C1_path).get_data()
            COPE[:,:,:,s] = C1
            print(subject)
        # save as nifti file with correct header and data type        
        outData = nib.Nifti1Image(COPE, affine=mni.affine, header=mni.header) # pass affine and header from last MNI image
        outData.set_data_dtype(np.float32)
        nib.save(outData, outFile)
        
        # open higher level job and write command as new line
        randomise_output = os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_cope1'.format(task))
        cmd = 'randomise -i {} -o {} -m {} -1 -x -c 3.1'.format(outFile,randomise_output,mask) #outFile is from concatenation
        self.higher_level_job = open(self.higher_level_job_path, "a") # append is important, not write
        self.higher_level_job.write(cmd)   # feat command
        self.higher_level_job.write("\n\n")  # new line
        self.higher_level_job.close()
        print('success: localizers_randomise_input')
        
    def rsa_letters_ev_conditions(self,task='rsa'):
        # outputs a dataframe with the letter and color conditions (general) for the EVs in the rsa-letters analysis
        # a-z black, then a-z color, 53rd EV was oddballs (nuisance regressor)

        #############################
        # output ev-number to letter-color condition file
        # make letter and color_condition columns based on EV numbers
        LCDF = pd.DataFrame()
        LCDF['ev'] = np.arange(1,53)
        LCDF['letter'] = self.letters*2        
        LCDF['trial_type_color'] = np.concatenate([np.repeat('black',len(self.letters)),np.repeat('color',len(self.letters))])
        LCDF.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_letters_ev_conditions.tsv'.format(task)),sep='\t')
        print('success: rsa_letters_ev_conditions')
        
    def rsa_letters_conditions(self,task='rsa'):      
        # concatenates all subjects colors files into a single one
        
        # load evs conditions file
        LCDF = pd.read_csv(os.path.join(self.higher_level_dir,'task-{}'.format(task),'task-{}_letters_ev_conditions.tsv'.format(task)),sep='\t')
        try:
            LCDF = LCDF.loc[:, ~LCDF.columns.str.contains('^Unnamed')]
        except:
            pass
        
        DFOUT = pd.DataFrame()
        for s,subject in enumerate(self.subjects):
            events = pd.read_csv(os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}_{}_events.tsv'.format(task,subject,'ses-01')),sep='\t')
            try:
                events = events.loc[:, ~events.columns.str.contains('^Unnamed')]
            except:
                pass
            # drop trial-wise information
            events = events[events['oddball']==0]
            colors = events.get(['subject','letter', 'trial_type_letter', 'trial_type_color','r','g', 'b',])
            # drop duplicates
            colors = colors.drop_duplicates()
            # a-z black, then a-z color
            colors = colors.sort_values(by=['trial_type_color','letter']).reset_index()
            # merge with ev condition file
            colors = pd.merge(colors, LCDF, on=["letter", "trial_type_color"])
            # concate
            DFOUT = pd.concat([DFOUT,colors],axis=0)
        # shell()
        DFOUT = DFOUT.loc[:, ~DFOUT.columns.str.contains('index')]
        DFOUT.to_csv(os.path.join(self.higher_level_dir,'task-{}'.format('rsa'),'task-rsa_letters_conditions.tsv'), sep='\t')
        print('success: rsa_letters_combine_colors')
    

        
    def housekeeping_dcm(self,task='rsa'):
        # For each subject and session of the RSA task, unzips then splits the nifti files into the DCM folder
        # copies the events file into the DCM folder
        # fslsplit <input> [output_basename] [-t/x/y/z]
        # gunzip all files in directory
        # input = first_level/task-rsa
        # output = AI project folder/Faye/DCM
        
        DCM_path = '/project/2422045.01/faye/DCM'
        dcm_sessions = ['session1','session2'] # different folder names
        
        # write unix commands to job to run in parallel
        self.dcm_housekeeping_path = os.path.join(self.analysis_dir,'jobs','job_dcm_fslsplit.txt')
        self.dcm_housekeeping = open(self.dcm_housekeeping_path, "w")
        self.dcm_housekeeping.write("#!/bin/bash\n")
        self.dcm_housekeeping.close()
        
        for s,subject in enumerate(self.subjects):
            for sess,session in enumerate(self.sessions):
                
                # make subject's directory in DCM folder path
                subj_dir_path = os.path.join(DCM_path, dcm_sessions[sess],'sub-{}'.format(subject))
                if not os.path.isdir(subj_dir_path):
                    os.mkdir(subj_dir_path)
                    
                #### FSL SPLIT ####
                # split it in the time dimension
                input_split = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}_{}_bold_mni.nii.gz'.format(task,subject,session))
                split_base = os.path.join(DCM_path, dcm_sessions[sess],'sub-{}'.format(subject),'task-{}_sub-{}_{}_bold_mni_'.format(task,subject,session))
                cmd1 = 'fslsplit {} {} -t'.format(input_split, split_base)
                
                #### UNZIP ####
                # first unzip all files in the subject's folder
                cmd2 = 'gunzip {}/*.gz'.format(subj_dir_path)
                
                #### copy events file into DCM folder ####
                old_events = os.path.join(self.first_level_dir,'task-{}'.format(task),'task-{}_sub-{}_{}_events.tsv'.format(task,subject,session))
                new_events = os.path.join(DCM_path, dcm_sessions[sess],'sub-{}'.format(subject),'task-{}_sub-{}_{}_events.tsv'.format(task,subject,session))
                cmd3 = 'cp {} {}'.format(old_events,new_events)
                
                # open job and write command as new line
                self.dcm_housekeeping = open(self.dcm_housekeeping_path, "a") # append is important, not write
                self.dcm_housekeeping.write(cmd1)   # split command
                self.dcm_housekeeping.write("\n\n")  # new line
                self.dcm_housekeeping.write(cmd2)   # gunzip command
                self.dcm_housekeeping.write("\n\n")  # new line
                self.dcm_housekeeping.write(cmd3)   # copy command
                self.dcm_housekeeping.write("\n\n")  # new line
                self.dcm_housekeeping.close()

        print('success: housekeeping_dcm')
        
    def qualia(self,):    
        # calculate the PA score from the Reading Questionnaire
        # see PAscoring 2020.xls for details
        # PA score = projector - associator
        
        flip_scores = ['7', '15', '5', '6', '8', '16'] # flip the response direction of the Likert scale from 5 to 1 for these questions
        # associator = ['7', '15', '5', '6', '8', '16', '10', '3', '14', '12', '18', '1']
        # projector = ['7', '15', '5', '6', '8', '16', '17', '2', '9', '11', '4', '13']
        
        associator = [ '10', '3', '14', '12', '18', '1']
        projector = ['17', '2', '9', '11', '4', '13']
                
        DF = pd.read_csv(os.path.join(self.higher_level_dir,'participants_qualia.csv'))
        
        DF2 = DF.copy()
        DF2[flip_scores] = np.max(DF[flip_scores]) + 1 - DF[flip_scores] # flip the likert scale
        
        DF2['associator'] = np.mean(DF[associator],axis=1)
        DF2['projector'] = np.mean(DF[projector],axis=1)
        DF['PA_score'] = DF2['projector'] - DF2['associator']
        DF.to_csv(os.path.join(self.higher_level_dir,'participants_qualia.csv'))
        # PA = projector - associator
    

