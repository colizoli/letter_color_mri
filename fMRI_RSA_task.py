#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
RSA Letter-Color Training, MRI
O.Colizoli & M.Blasco Oliver, 2019
Event related design - 2x2 - Trained vs. Untrained, Color vs. Black
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
REPS       = 2        # how many times to repeat all trials in trials_RSA
        
## load CSV file for balancing stimuli
trials_RSA = pd.read_csv(os.path.join('stimuli','trials_RSA.csv'),sep='\t',header=None)
#print(trials_RSA)
trials_RSA.columns=['letter_idx','oddball','letter_condition','color_condition']
trials_RSA.loc[0,'letter_idx']=0 # this version of psychopy adding a "?" to first row always, but not to S1.csv file...


# Odd stimuli
oddball_numbers = np.array(['1','2','3','4','5','6','7','8','9'])
oddball_colors = [
    [128,128,128],
    [188,188,188],
    [117,117,117]
]

# Timing (seconds):
t_stim = 1.5  # stimulus presentation 
t_ITI  = [4,6]   # jitter, uniform
t_RT   = 3

#########################################
################# TEXTS #################
#########################################
welcome_txt = "You will see a series of letters and numbers, presented in black, grey or color.\
\n\nYour task is to respond to 'targets':\
\ni) ANY of these NUMBERS: 1,2,3,4,5,6,7,8,9\
\nii) a GREY COLOR (letter or number!)\
\nWhenever you see either i) or ii),\
\nthen you must press the button with your RIGHT INDEX FINGER to respond.\
\n\nAt the end of the scan, you will be told how well you performed.\
\n\nMaintain eye-fixation at the center of the screen at all times (even when there are no letters or numbers).\
\nLay as still as possible for the duration of the scan session (that includes between scans!)\
\n\n[Waiting for scanner...]"

#########################################
############## MAIN FUNCS ###############
#########################################

# Get subject number
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
# run = 1

if subject_ID:
    
    ## LETTER SET DEPENDENT ON SUBJECT NUMBER
    untrained_letters = gp.letter_sets_lower[not np.mod(subject_ID,2)] # lower case, UNTRAINED set!!
    
    # Get subject-specific letter-color mapping, letters always have same order, colors are permuted
    DFS = pd.read_csv(os.path.join('stimuli','S{}.csv'.format(subject_ID)))
    print(DFS)
    subj_colors = [[DFS['r'][i],DFS['g'][i],DFS['b'][i]] for i in range(len(DFS))]
    subj_letters = list(DFS['letter'])
    
    # concat, trained + untrained, then zip
    letters = subj_letters + untrained_letters
    colors = subj_colors*2
    mappings = list(zip(letters,colors))       #1-1 mapping
        
    ## Create LogFile folder cwd/LogFiles
    cwd = os.getcwd()
    logfile_dir = os.path.join(cwd,'LogFiles')  # system independent (Marta's way doesn't work on Mac)
    
    if not os.path.isdir(logfile_dir):
        os.makedirs(logfile_dir)
        
    ## output file name with time stamp prevents any overwriting of data
    timestr = time.strftime("%Y%m%d-%H%M%S") 
    output_filename = os.path.join(logfile_dir,'RSA_S{}_sess{}_run{}_{}.csv'.format(subject_ID,session,run,timestr))
    cols = [
        'subject',
        'session',
        'run',
        'trial_num',
        'letter',
        'letter_condition',
        'color_condition',
        'r','g','b',
        'ITI',
        'oddball',
        'button_pressed',
        'correct',
        'RT',
        'onset_time'
    ]

    DF = pd.DataFrame(columns=cols) 

    # Set-up window:
    mon = monitors.Monitor('myMac15', width=p.screen_width, distance=p.screen_dist)
    mon.setSizePix((p.scnWidth, p.scnHeight))
    win = visual.Window((p.scnWidth, p.scnHeight),color=gp.white,colorSpace='rgb255',monitor=mon,fullscr=not debug_mode,units='pix',allowStencil=True,autoLog=False)
    win.setMouseVisible(False)

    # Set-up stimuli and timing
    instr         = visual.TextStim(win, color='black', pos=(0.0, 0.0), wrapWidth=p.ww ) 
    stim_fix      = visual.TextStim(win, text='+', color='black', pos=(0.0, 0.0), height=p.fh) 
    stim_letter   = visual.TextStim(win, color ='black', pos=(0.0, p.sp), height=p.lh2, font=gp.font_trained) 
    clock         = core.Clock()
    
    # Set conditions and stimulus list    
    ALL_TRIALS  = pd.concat([trials_RSA]*REPS, ignore_index=True)    
    ALL_TRIALS = ALL_TRIALS.sample(frac=1).reset_index(drop=True) ## SHUFFLE ORDER OF TRIALS
    if debug_mode:
        ALL_TRIALS = ALL_TRIALS.iloc[0:13,:] # just get first X trials
    # print(ALL_TRIALS)
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
    
    correct_resp = [] # for feedback at 
    targets_resp = [] 
    #### MAIN LOOP ###
    for t in range(len(ALL_TRIALS)):

        # First, check if oddball or not
        this_answer = int(ALL_TRIALS['oddball'][t])
        this_letter_condition = ALL_TRIALS['letter_condition'][t]
        this_color_condition = ALL_TRIALS['color_condition'][t]
        if this_answer: # ODDBALL
            # Check letter or number?
            if int(ALL_TRIALS['letter_idx'][t]) > 26: # NUMBER
                this_letter = oddball_numbers[random.randint(0, 8)]  # choose random 1-9
            else:
                this_letter = letters[random.randint(0,25)] # choose random letter from whole alphabet
            # Check color or grey?    
            if  this_color_condition == 'color': # COLOR
                this_color =  colors[random.randint(0,25)]  
            else:   # GREY
                this_color =  oddball_colors[random.randint(0,2)]  # greys to choose from           
        else: # NOT oddball
            this_letter = mappings[int(ALL_TRIALS['letter_idx'][t])][0]
            # Check color condition
            if this_color_condition == 'color': # COLOR
                this_color =  mappings[int(ALL_TRIALS['letter_idx'][t])][1]
            else:
                this_color = gp.black
            
        # Target stimuli current trial        
        print('########## Trial {} #########'.format(t+1))
        print('Color of current trial: {}'.format(this_color))
        print('Letter of current trial: {}'.format(this_letter))
        print(this_answer)
        
        # STIMULUS (letter or number, black or color) & wait for response if odd ball
        stim_letter.setText(this_letter)
        stim_letter.setColor(this_color,'rgb255')
        stim_letter.draw()
        win.flip()
        stim_onset = clock.getTime()
        core.wait(t_stim)
        # Inter-trial interval, jittered
        ITI = np.round(random.uniform(t_ITI[0],t_ITI[1]),2) # draw from uniform distribution
        stim_fix.draw()
        win.flip()
        core.wait(ITI)
        
        response = event.getKeys(keyList=p.buttons+['q'],timeStamped=clock) # MONITOR FOR POSSIBLE KEY PRESS
        if response:
            button_pressed = response[0][0]
            # q quits the experiment
            if button_pressed == 'q':
                core.quit()
            else:
                RT = response[0][1] - stim_onset # latency relative to clock.reset()
                correct = (response[0][0] in p.buttons) and this_answer # odd ball and pressed button
                targets_resp.append(correct)
        else: 
            button_pressed = np.nan   # didn't press
            correct = not this_answer # not an oddball
            RT = np.nan
        correct_resp.append(correct) # save correct responses
        
        # output data frame on each trial 
        DF.loc[t] = [
            int(subject_ID),    # subject
            session,            # session
            run,                # run
            t,                  # trial_number
            this_letter,        # frequency
            this_letter_condition,     # letter condition
            this_color_condition,       # color condition
            int(this_color[0]), # r
            int(this_color[1]), # g
            int(this_color[2]), # b
            ITI,                # ITI
            int(this_answer),   # correct answer, 1=oddball,0=normal
            button_pressed,     # button pressed
            int(correct),       # correct response or error?
            round(RT,8),        # RT
            round(stim_onset-first_pulse,8) # onset_time relative to first pulse
        ]
        DF.to_csv(output_filename)
                     
    # BOLD BASELINE 
    stim_fix.draw()
    win.flip()
    core.wait(p.t_baseline)

# Scan time
print('Scan time (s) = {}'.format(clock.getTime() - first_pulse))

## PRESENT ACCURACY ## 
avg_acc = np.true_divide(np.sum(correct_resp),len(ALL_TRIALS))*100
ntargets = np.sum(ALL_TRIALS['oddball'])
nhits = np.sum(targets_resp)
txt = 'Done! Your total accuracy on this block was {}% correct.\
\nYou hit {} out of {} targets'.format(round(avg_acc,2), nhits, ntargets)
instr.setText(txt)
instr.draw()
win.flip()
core.wait(5) 

# Close-up 
win.close()
core.quit()


