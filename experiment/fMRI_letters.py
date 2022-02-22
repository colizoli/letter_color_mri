#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
Single letter localizer, MRI
O.Colizoli & M.Blasco Oliver, 2019
ABAB blocked design
Short rest period between each block (fixation)
BOLD baseline period beginning and end of the run (fixation)
"""

# Import necessary modules
import os, time # for paths and data
from psychopy import core, visual, event, gui, monitors
import random
import numpy as np
import pandas as pd
import serial
#from IPython import embed as shell
import mri_parameters as p # see for timing, stimuli and parameters
import general_parameters as gp # for letter sets, counterbalancing

#########################################
##############   PARAMS   ###############
#########################################
debug_mode = False

#########################################
################# TEXTS #################
#########################################
welcome_txt = "You will see a series of letters and symbols presented on screen.\
\nYour task is to (passively) view the letters and symbols,\
\nyou do not have to respond.\
\n\nMaintain eye-fixation at the center of the screen at all times\
\n(even when there are no letters or symbols).\
\n\nLay as still as possible for the duration of the scan session\
\n(that includes between the scans!)\
\n\n[Waiting for scanner...]"

#########################################
############## MAIN FUNCS ###############
#########################################

#Get subject number
g = gui.Dlg()
g.addField('Subject Number:')
g.addField('Session:')
g.addField('Run:')
g.show()
subject_ID = int(g.data[0])
session = int(g.data[1])
run = int(g.data[2])
# in case GUI doesn't work
# subject_ID = 1
# session = 1

if subject_ID:
    
    ## LETTER SET DEPENDENT ON SUBJECT NUMBER (UNTRAINED SET)
    letters = gp.letter_sets_upper[not np.mod(subject_ID,2)] # upper case, UNTRAINED set!!
    
    ## Create LogFile folder cwd/LogFiles
    cwd = os.getcwd()
    logfile_dir = os.path.join(cwd,'LogFiles','sub-{}'.format(subject_ID),'sess-{}'.format(session),'func')  # system independent
    
    if not os.path.isdir(logfile_dir): # make the folder if it doesn't exist already
        os.makedirs(logfile_dir)
        
    ## output file name with time stamp prevents any overwriting of data
    timestr = time.strftime("%Y%m%d-%H%M%S") 
    # bids format    
    output_filename = os.path.join(logfile_dir,'sub-{}_sess-{}_task-lettersloc_run-{}_events_{}.tsv'.format(subject_ID,session,run,timestr))
    cols = ['subject',
              'session',
              'block',
              'trial_type',
              'trial_number',
              'font',
              'letter',
              'onset',
              'duration']
      
    DF = pd.DataFrame(columns=cols) #called DF in other script  

    # Set-up window:
    mon = monitors.Monitor('myMac15', width=p.screen_width, distance=p.screen_dist)
    mon.setSizePix((p.scnWidth, p.scnHeight))
    win = visual.Window((p.scnWidth, p.scnHeight),color=gp.white,colorSpace='rgb255',monitor=mon,fullscr=not debug_mode,units='pix',allowStencil=True,autoLog=False)
    win.setMouseVisible(False)

    # Set-up stimuli and timing
    instr       = visual.TextStim(win, color='black', pos=(0.0, 0.0), wrapWidth=p.ww ) 
    stim_fix    = visual.TextStim(win, text='+', color='black', pos=(0.0, 0.0), height=p.fh) 
    stim_word   = visual.TextStim(win, color = 'black', pos=(0.0, 0.0), height=p.lh) # SET FONT HERE, SCALING
    stim_word.fontFiles = [os.path.join('stimuli','HebrewUniversal.ttf')]
    stim_word.fontFiles = [os.path.join('stimuli','Tahoma.ttf')]
    clock       = core.Clock()
    
    # Set conditions and stimulus list
    if np.mod(subject_ID,2): # counterbalance letter-symbol order
        p.blocks_letter.reverse()
    blocks          = p.blocks_letter*p.reps_loc # total number of blocks
    trials          = p.trials_loc              # per block
    if debug_mode:
        blocks      = p.blocks_letter*1
        
    ##################################
    # Instructions, waiting for scanner
    clock.reset() # important NOT to reset clock after this!
    instr.setText(welcome_txt)
    instr.draw()
    win.flip()
    # WAIT FOR TRIGGER
    # respond = event.waitKeys(keyList=TRIGGER_CODE, timeStamped=clock)
    respond = p.waitForTrigger(clock)
    key_pressed, first_pulse = respond[0] # in seconds
    ##################################

    # BOLD BASELINE 
    stim_fix.draw()
    win.flip()
    # For quitting early
    respond = event.waitKeys(maxWait=p.t_baseline, keyList=['q'])
    if respond:
        key_pressed, latency = respond[0]
        if key_pressed == 'q':
            #win.close() #added my Marta, otherwise windows crash
            core.quit()

    counter_main = 0 
    #### MAIN LOOP ###
    for b in range(len(blocks)):
        condition   = blocks[b]
        # shuffle letters
        letters    = pd.DataFrame(letters)
        letters    = letters.sample(frac=1).reset_index(drop=True) # SHUFFLE LIST
        
        # trials per block   
        for t in range(trials):        
            # Define stimuli and variables according to block type
            this_word          = letters.loc[t][0]
            if condition == 'Letter':     # letter
                this_font      = p.fonts_letter
                this_size      = 1                  # scaling of fonts
                this_pos       = 0.0                # y position of font
                
            elif condition == 'Symbol':   # symbol
                this_font      = p.fonts_symbol
                this_size      = p.ss
                this_pos       = p.sp               # y position of fonts
                
            #Show fixation cross
            stim_fix.draw() 
            win.flip()
            core.wait(p.t_fix)

            #Present word
            stim_word.setText(this_word)
            stim_word.setFont(this_font)
            stim_word.setHeight(p.lh*this_size) # scaling of fonts
            stim_word.setPos((0.0,this_pos))  # y position of fonts
            stim_word.draw() 
            win.flip()
            this_word_time = clock.getTime() - first_pulse # get timing of word with respect to first pulse
            core.wait(p.t_word)

            # Save trial data       
            DF.loc[counter_main] = [subject_ID,     # Subject number
                                session,            # Session
                                b+1,                # Block number
                                condition,          # Language
                                counter_main+1,     # trial number in the whole task
                                this_font,          # Font used in current trial
                                this_word,          # Word presented in current trial
                                this_word_time,     # Word onset with respect to first_pulse
                                p.t_word        ]   # duration (seconds)  
                                    
            DF.to_csv(output_filename,sep='\t')
            counter_main = counter_main+1
            
            # For quitting early
            keys = event.getKeys()
            if keys:
                # q quits the experiment
                if keys[0] == 'q':
                    core.quit()
                    
        # Time between blocks
        stim_fix.draw() 
        win.flip()
        core.wait(p.IBI)
          
    # BOLD BASELINE 
    stim_fix.draw()
    win.flip()
    core.wait(p.t_baseline)

# Scan time
print('Scan time (s) = {}'.format(clock.getTime() - first_pulse))

# Close-up 
stim_fix.setColor(gp.white,'rgb255')  
stim_fix.draw()
win.flip()
core.wait(1.5) # just wait max 1 TR for scanner to stop
win.close()
core.quit()














