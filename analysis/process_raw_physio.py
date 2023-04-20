#!/usr/bin/env python
# encoding: utf-8

"""
process_raw_physio.py
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
The "position in data points" (idx) column in the .vmrk file with the trigger markers
corresponds to the index (row) of the data file (.eeg). The "time" column in df is 
irrelevant.

In the .vhdr files, the following have to match the filenames themselves:
e.g.,
DataFile=sub-211_sess-1.eeg
MarkerFile=sub-211_sess-1.vmrk

change:
sub-211_sess-1
sub-201_sess-1
sub-136_sess-1
sub-103_sess-1
sub-109_sess-1
"""

import os
import glob
import mne
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
SAMPLE_ERROR = 50 # margin of error in markers (Hz)
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
    """Iterate through the run_gap and split the dataframe into multiple .csv files
    
    Args:
        filename (str): the path to the .vhdr file
        df (pandas.core.frame.DataFrame): the physio data dataframe (.eeg)
        df_triggers (pandas.core.frame.DataFrame): the trigger marker dataframe
    
    Notes:
    ------
        pandas.core.frame.DataFrames with index, heart rate, and respiration saved (cols) for each run
        A run will get all data until the next run starts (i.e., there will be extra data on tail ends)
    """
    
    df.drop(['time'],axis=1,inplace=True) # time irrelevant, only row position important
    
    # The run starts are where there are large run_gaps and first trigger
    df_triggers['starts'] = np.array((df_triggers['run_gap'] > int((SAMPLE_RATE * TR) + SAMPLE_ERROR)) | (df_triggers['run_gap'] == 0), dtype=int)
    
    # get_stops = np.array(np.diff(np.array(df_triggers['starts'])) > 0, dtype=int)
    # df_triggers['stops'] =  np.append(get_stops, [0])
    
    df_starts = df_triggers[df_triggers['starts']==1].copy() # get only starts 
    
    # loop through starting indexes, get data out of df
    start_idx = list(df_starts['idx']) # list of starting indexes
    for run_counter,idx in enumerate(start_idx):
        if idx == start_idx[-1]: # check if last run
            # go till end of dataframe
            this_run = df.iloc[int(idx):,:]            
        else: # otherwise, get data until next starting index
            this_run = df.iloc[int(idx):int(start_idx[run_counter+1]),:]
            
        # save run in format to sub-000_ses-01_run-01_physio)
        fn_out = os.path.join(output_dir, os.path.splitext(filename)[0] + '_run-{:02}_physio.tsv'.format(run_counter+1))
        
        # housekeeping
        try:
            fn_out = fn_out.replace('sess-1', 'ses-01')
            fn_out = fn_out.replace('sess-2', 'ses-02')
        except:
            pass
            
        this_run.to_csv(fn_out, sep='\t') 
        print('success: split_dataframe {}'.format(fn_out))
        
        
def loop_raw_physio():
    """Loop through all files in "physiology" and run workflow.
    """
    for fn in os.listdir(raw_dir):
        if '.vhdr' in fn:
            df = convert_data(os.path.join(raw_dir,fn)) # convert to MNE dataframe
            df_triggers = find_triggers(os.path.join(raw_dir,fn)) # get trigger markers dataframe
            split_dataframe(fn, df, df_triggers) # split data into runs


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
                
# -----------------------
# Run
# -----------------------   
if __name__ == "__main__":
    # loop_raw_physio()
    descriptives_physio()
