#!/usr/bin/env python
# encoding: utf-8

"""
process_source_physio.py
Created by O.Colizoli, C. de Witte April-2023
Last update: 20-04-2023
Python version 3.6
************************
MRI letter-color project
************************
A single time series of heart rate and respiration data was recorded 
during scanning using BrainVision Recorder. 
This code splits the time series into the runs based on the MRI triggers. 
Once the runs are split, they will have to be matched to the scan 
acquisition to remove unfinished/uncompleted/unneeded runs.


Notes:
------
In the final text files, the triggers have to have a separate column!
Save final text finals with extension .txt
Columns of output: 'idx', 'HR', 'Resp', 'trigger' but no column headers 

The "position in data points" (idx) column in the .vmrk file with the trigger markers
corresponds to the index (row) of the data file (.eeg). The "time" column in df is 
irrelevant.

In the .vhdr files, the following have to match the filenames themselves:
e.g.,
DataFile=sub-211_sess-1.eeg
MarkerFile=sub-211_sess-1.vmrk

change:
sub-103_sess-1
sub-109_sess-1
sub-119_sess-1
sub-136_sess-1
sub-201_sess-1
sub-211_sess-1
sub-213_sess-2
"""

import os
import glob
import mne
import shutil as sh
import pandas as pd
import numpy as np
from IPython import embed as shell # for Oly's debugging only

# -----------------------
# Paths
# ----------------------- 
##########################
# home_dir = os.path.dirname(os.getcwd()) # one level up from analysis folder
home_dir = '/project/3018051.01/'
##########################
raw_dir    = os.path.join(home_dir, 'bids_raw','physiology')  # raw physiology files
output_dir = os.path.join(home_dir, 'bids_raw','physiology_processed') # after split
subject_dir    = os.path.join(home_dir, 'analysis')  # subject list
bids_dir       = os.path.join(home_dir, 'derivatives')        # NIFTI versions of DICOM files for processing

if not os.path.isdir(output_dir): # make the output directory if doesn't exist
    os.mkdir(output_dir)

# -----------------------
# Scan trigger in marker file (.vmrk)
# -----------------------     
MARKER_TRIGGER = 'R  1'
# -----------------------
# Sampling rate (Hz)
# -----------------------     
SAMPLE_RATE = 5000
SAMPLE_ERROR = 5 # margin of error in markers (Hz)
# -----------------------
# Repetition Time (s)
# -----------------------   
TR = 1.5

# -----------------------
# Define functionality
# -----------------------     
def convert_data(filename):
    """Read header file with MNE-python and turn into dataframe
    
    Args:
        filename (str): the path to the .vhdr file
    
    Returns:
        pandas.core.frame.DataFrame: raw data converted to pandas datafrme
    """
    raw = mne.io.read_raw_brainvision(filename, preload=True)
    df = raw.to_data_frame()
    
    return df
    
    print('success: convert_data {}'.format(filename))


def find_triggers(filename):
    """Read the marker text file and return trigger markers.
    
    Args:
        filename (str): the path to the .vhdr file
    
    Returns:
        pandas.core.frame.DataFrame: scan trigger markers converted to pandas datafrme
    """
    with open(os.path.splitext(filename)[0] + '.vmrk', 'r') as f:
        for i in range(11): # skip first 11 lines in header file
            next(f)
        markers = [line.strip().split(',') for line in f]
    
    # Turn text into a dataframe
    mark_df = pd.DataFrame(markers, columns=['mk', 'description', 'idx', 'size', 'chnr', 'col6'])
    
    # Find the rows where description equals the marker_trigger
    mask = mark_df['description'] == MARKER_TRIGGER
    
    # get only those rows with MARKER_TRIGGER
    df_triggers = mark_df[mask].copy()
    
    # drop useless columns
    df_triggers.drop(['size', 'chnr', 'col6'], axis=1, inplace=True)
    
    # check if distance between triggers is greater than one repetition time (TR)
    run_gap = np.diff(np.array(df_triggers['idx'],dtype=int))
    
    # add 0 to the first element to make same length as df_triggers again
    df_triggers['run_gap'] = np.insert(run_gap, 0, [0])

    return df_triggers # Return the resulting dataframe with triggers and run_gap


def split_dataframe(filename, df, df_triggers):
    """Iterate through the run_gap and split the dataframe into multiple .csv files.
    
    Args:
        filename (str): the path to the .vhdr file
        df (pandas.core.frame.DataFrame): the physio data dataframe (.eeg)
        df_triggers (pandas.core.frame.DataFrame): the trigger marker dataframe
    
    Notes:
    ------
        pandas.core.frame.DataFrames with index, heart rate, and respiration saved (cols) for each run
        A run will get data from start to stop+TR (at least 3 TRs distance, because sometimes 1 trigger is missed)
    """
    
    df.drop(['time'],axis=1,inplace=True) # time irrelevant, only row position important
    df_triggers['idx'] = pd.to_numeric(df_triggers['idx'])
    
    # The run starts are where there are large run_gaps (more than 3 TRs) and first trigger
    df_triggers['starts'] = np.array((df_triggers['run_gap'] >= int(SAMPLE_RATE * TR * 3)) | (df_triggers['run_gap'] == 0) , dtype=int)
        
    # get stops
    get_stops = np.array(np.diff(np.array(df_triggers['starts'])) > 0, dtype=int)
    df_triggers['stops'] =  np.append(get_stops, [1]) # last trigger is always a stop
    
    df_starts = df_triggers[df_triggers['starts']==1].copy() # get only starts 
    df_stops = df_triggers[df_triggers['stops']==1].copy() # get only stops
    
    start_idx = list(df_starts['idx']) # list of starting indexes
    stop_idx = list(df_stops['idx']) # list of stopping indexes
    
    # loop through starting indexes, get data out of df and triggers out of df_triggers
    for run_counter,idx in enumerate(start_idx):
        
        idx_start = int(idx)
        idx_stop  = int(stop_idx[run_counter])+int(SAMPLE_RATE * TR)
        
        # get data until stop + TR 
        this_run = df.iloc[idx_start:idx_stop,:]        
        
        # select triggers only within the min-max indices of subj_funs...
        this_triggers = df_triggers[df_triggers['idx'].between(idx_start, idx_stop)]
        this_triggers = pd.DataFrame(this_triggers['idx'])
        this_triggers['trigger'] = np.repeat(1,this_triggers.shape[0])
                
        this_run.reset_index(inplace=True)
        this_run.columns = ['idx', 'HR', 'Resp']        
        
        subj_merged = this_run.merge(this_triggers, on='idx', how='outer', indicator=False) # outer merge
                
        subj_merged['trigger'] = subj_merged['trigger'].fillna(0) # nan to 0
        subj_merged['trigger'] = subj_merged['trigger'].astype('int')        

        # save run in format to sub-000_ses-00_run-00_physio.tsv
        path_subj = os.path.splitext(filename)[0] 
        # housekeeping
        try:
            path_subj = path_subj.replace('sess-1', 'ses-01')
            path_subj = path_subj.replace('sess-2', 'ses-02')
        except:
            pass
        fn_out = os.path.join(output_dir, path_subj + '_run-{:02}_physio.tsv'.format(run_counter+1))
        
        # drop index and column headers
        subj_merged.to_csv(fn_out, sep='\t', header=False, index=False) 
        
        print('splitting... {}'.format(fn_out))
    # save triggers too
    df_triggers.to_csv(os.path.join(output_dir, path_subj + '_triggers.csv'))
    print('success: split_datarame')
    
            
def loop_raw_physio():
    """Loop through all files in "physiology" and run workflow.
    """
    for fn in os.listdir(raw_dir):
        if '.vhdr' in fn:
            df = convert_data(os.path.join(raw_dir,fn)) # convert to MNE dataframe
            df_triggers = find_triggers(os.path.join(raw_dir,fn)) # get trigger markers dataframe
            split_dataframe(fn, df, df_triggers) # split data into runs
    print('success: loop_raw_physio')


def descriptives_physio():
    """Count the number of runs and TRs for each subject and session.
    
    Notes:
    ------
        pandas.core.frame.DataFrames output with subject, session, runs, samples and TRs.
    """
    
    df_out = pd.DataFrame(columns=['subject','session','runs','samples','trs'])
    
    dfs = pd.read_csv(os.path.join(subject_dir, 'participants_full_mri.csv'))

    counter = 0
    for s, subj in enumerate(dfs['subjects']):
        for sess in ['ses-01', 'ses-02']:
            
            subj_fns = glob.glob(os.path.join(output_dir, 'sub-{}_{}_run-*_physio.tsv'.format(subj,sess)))
            
            # loop through all runs for current subject, session
            for r,fn in enumerate(subj_fns):    

                df_phys = pd.read_csv(fn, sep='\t')
                
                # output row of data frame for each run 
                df_out.loc[counter] = [
                    int(subj),                          # subject
                    sess,                               # session
                    int(r+1),                           # run
                    int(len(df_phys)),                  # no. samples
                    int(len(df_phys)/(SAMPLE_RATE*TR)), # no. TRs
                ]
                df_out.to_csv(os.path.join(output_dir, 'descriptives_physio.csv'))
                counter += 1
            print('descriptives_physio: subject {} session {}'.format(subj, sess))
    print('success: descriptives_physio')


def rename_physio():
    """Rename extracted physiological runs to psychopy tasks.
    
    Notes:
    ------
        The localizers had 248 TRs, the RSA runs had 496 or 497 TRs
        Match from end of session: LOC2 first, RSA1 last
        When ready copy the files from the main script: run_letter_color_analysis.py
    """
    df_all = pd.read_csv(os.path.join(output_dir, 'descriptives_physio.csv'))
    
    runs = ['task-loc2', 'task-loc1', 'task-rsa_run-04', 'task-rsa_run-03', 'task-rsa_run-02', 'task-rsa_run-01']        
    
    # loop over subjects and copy to correct functional folder
    dfs = pd.read_csv(os.path.join(subject_dir, 'participants_full_mri.csv'))
    
    flag_subjects = [] # to correct by hand
    
    for s, subj in enumerate(dfs['subjects']):
        for sess in ['ses-01', 'ses-02']:
            subj_fns = glob.glob(os.path.join(output_dir, 'sub-{}_{}_run-*_physio.tsv'.format(subj,sess)))
            
            if not subj == 201: # localizer first

                subj_fns = subj_fns[::-1] # started with localizer
                if len(subj_fns) == 6:
                    for rcounter, phys in enumerate(subj_fns):
                        orig_fname = phys.split("/")[-1].split("_")

                        old_name = phys
                        new_name = os.path.join(output_dir, '{}_{}_{}_physio.tsv'.format(orig_fname[0], orig_fname[1], runs[rcounter]))
                        print(old_name)
                        print(new_name)
                        print()
                        
                        os.rename(old_name, new_name)
                        # cmd = sh.copyfile(src, dst)
                else:
                    #flag
                    flag_subjects.append('sub-{} {}'.format(subj, sess))
    print('WARNING: check the following by hand and correct')
    print(flag_subjects)
    print('success: rename_physio')

    
# -----------------------
# Run
# -----------------------   
if __name__ == "__main__":
    loop_raw_physio()
    descriptives_physio()
    rename_physio()

    