#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
Color localizer, MRI
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
from IPython import embed as shell
import mri_parameters as p # see for timing, stimuli and parameters
import general_parameters as gp # for letter sets, counterbalancing

#########################################
##############   PARAMS   ###############
#########################################
debug_mode = False

#########################################
################# TEXTS #################
#########################################
welcome_txt = "You will see a series of symbols presented on screen.\
\nYour task is to (passively) view the symbols,\
\nyou do not have to respond.\
\n\nMaintain eye-fixation at the center of the screen at all times\
\n(even when there are no symbols).\
\n\nLay as still as possible for the duration of the scan session\
(that includes between scans!)\
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
    
    ## Create LogFile folder cwd/LogFiles
    cwd = os.getcwd()
    logfile_dir = os.path.join(cwd,'LogFiles','sub-{}'.format(subject_ID),'sess-{}'.format(session),'func')  # system independent
    
    if not os.path.isdir(logfile_dir): # make the folder if it doesn't exist already
        os.makedirs(logfile_dir)
        
    ## output file name with time stamp prevents any overwriting of data
    timestr = time.strftime("%Y%m%d-%H%M%S") 
    # bids format    
    output_filename = os.path.join(logfile_dir,'sub-{}_sess-{}_task-colorsloc_run-{}_events_{}.tsv'.format(subject_ID,session,run,timestr))
 
    cols = ['subject',
              'session',
              'block',
              'trial_type',
              'trial_number',
              'letter',
              'r',
              'g',
              'b',
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
    stim_word   = visual.TextStim(win, color ='black', pos=(0.0, p.sp), height=p.lh*p.ss) # SET FONT HERE, SCALING
    stim_word.fontFiles = [os.path.join('stimuli','HebrewUniversal.ttf')]
    stim_word.font = 'HebrewUniversal'
    clock       = core.Clock()
    
    
    # Set conditions and stimulus list
    if not np.mod(subject_ID,2): # counterbalance color-black order
        p.blocks_color.reverse()
    blocks          = p.blocks_color*p.reps_loc # total number of blocks
    trials          = p.trials_loc              # per block
    if debug_mode:
        blocks      = p.blocks_color*1
        
    ##################################
    # Instructions, waiting for scanner
    clock.reset() # important NOT to reset clock after this!
    instr.setText(welcome_txt)
    instr.draw()
    win.flip()
    # WAIT FOR TRIGGER
    # respond = event.waitKeys(keyList=p.trigger_key, timeStamped=clock)
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
            core.quit()

    counter_main = 0 
    #### MAIN LOOP ###
    for b in range(len(blocks)):
        condition   = blocks[b]
        # shuffle symbols and colors
        symbols    = pd.DataFrame(p.symbols)
        colors     = pd.DataFrame(p.colors_RGB)
        symbols    = symbols.sample(frac=1).reset_index(drop=True) # SHUFFLE LIST
        colors     = colors.sample(frac=1).reset_index(drop=True) # SHUFFLE LIST
        
        # trials per block   
        for t in range(trials):        
            # Define stimuli and variables according to block type
            this_word           = symbols.loc[t][0]
            if condition == 'Color':  # color
                this_color      = np.array(colors.loc[t])
            elif condition == 'Black':  # black
                this_color      = np.array([0,0,0]    )   
                            
            #Show fixation cross
            stim_fix.draw() 
            win.flip()
            core.wait(p.t_fix)

            #Present word
            stim_word.setText(this_word)
            stim_word.setColor(this_color,'rgb255')
            stim_word.draw() 
            win.flip()
            this_word_time = clock.getTime() - first_pulse # get timing of word with respect to first pulse
            core.wait(p.t_word)

            # Save trial data       
            DF.loc[counter_main] = [subject_ID,     # Subject number
                                session,            # Session
                                b,                  # Block number
                                condition,          # Language
                                counter_main,       # trial number in the whole task
                                this_word,          # Word presented in current trial
                                int(this_color[0]),      # R
                                int(this_color[1]),      # G
                                int(this_color[2]),      # B
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














