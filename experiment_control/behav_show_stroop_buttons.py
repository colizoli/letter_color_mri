"""
Train button-color responses for the behavioral stroop task
O.Colizoli 2023
"""
# data saved in ~/source/sub-XXX/ses-X

# Import necessary modules
import os, time  # for paths and data
import random
from psychopy import core, visual, event, data, sound, gui, monitors
import numpy as np
import pandas as pd
import general_parameters as p # for letter sets, counterbalancing
#from rusocsci import buttonbox
#from IPython import embed as shell

debug_mode = False

cwd = os.getcwd() # current working directory

stroop_letter_conditions = [
    ['e', 'n', 'o', 'r'], # counter balanced letters
    ['a', 'i', 's', 't' ]
]

## Keyboard Buttons 
buttons = ["q","q","q","q"] # colors 1,2,3,4 in order

## load CSV file for balancing stimuli
button_trials = [1,2,3,4]

# Get subject number
g = gui.Dlg()
g.addField('Subject Number:')
g.show()
subject_ID = int(g.data[0])
group = [3 if subject_ID < 200 else 4]

if subject_ID:
    ### Define stimuli for this subject & save file
    subject_colors  = os.path.join(cwd,'colors','sub-{}'.format(subject_ID),'sub-{}_colors.tsv'.format(subject_ID))
    DFS = pd.read_csv(subject_colors,sep='\t')
    # check which letter counterbalancing condition we have
    if np.sum(DFS['letter'].isin(stroop_letter_conditions[0])):
        DFS = DFS[DFS['letter'].isin(stroop_letter_conditions[0])] # select only high frequency letters for stroop task  
    else:
        DFS = DFS[DFS['letter'].isin(stroop_letter_conditions[1])] # select only high frequency letters for stroop task  
        
    DFS.sort_values('letter', ascending=True, inplace=True) # sort alphabetically to make sure the buttons are mapped correctly
    DFS.reset_index(inplace=True) # otherwise list comprehension fails
    subj_colors = [[DFS['r'][i],DFS['g'][i],DFS['b'][i]] for i in range(len(DFS))] # get RGB codes
    subj_letters = list(DFS['letter'])
    
    # Set-up window:
    mon = monitors.Monitor('myMac15', width=p.screen_width, distance=p.screen_dist)
    mon.setSizePix((p.scnWidth, p.scnHeight))
    win = visual.Window((p.scnWidth, p.scnHeight),color=p.white,colorSpace='rgb255',monitor=mon,fullscr=not debug_mode,units='pix',allowStencil=True,autoLog=False)
    win.setMouseVisible(False)
    
    # Set-up stimuli and timing
    welcome_txt = 'Learn Button-Color Associations! \
    \n\nYou will see colors appear one at a time on the screen.\
    \nInstructions: Press the correct button for each COLOR as fast AND as accurately as possible.\
    \nThe response options are shown in the colored squares below.\
    \nPlease use only one hand to respond - your dominant hand.\
    '

    stim_instr  = visual.TextStim(win, color='black', pos=(0.0, 0.0), wrapWidth=p.ww)  # can't really center, and TextBox doesn't work, stupid!
    stim_fix    = visual.TextStim(win, text='+', color='black', pos=(0.0, 0.0), height=p.fh)
    stim_sq     = visual.Rect(win, width=100, height=100, pos=(0.0, 0.0))
    feed_good   = visual.TextStim(win, text='Great!', color='green', pos=(0.0, 50.0))  # can't really center, and TextBox doesn't work, stupid!
    feed_error  = visual.TextStim(win, text='Wrong!', color='red', pos=(0.0, 50.0))  # can't really center, and TextBox doesn't work, stupid!
    feed_miss   = visual.TextStim(win, text='Too slow!', color='blue', pos=(0.0, 50.0))  # can't really center, and TextBox doesn't work, stupid!
    stim_sq1    = visual.Rect(win, width=50, height=50, pos=(-100.0, -200.0))
    stim_sq2    = visual.Rect(win, width=50, height=50, pos=(-50, -200.0))
    stim_sq3    = visual.Rect(win, width=50, height=50, pos=(50.0, -200.0))
    stim_sq4    = visual.Rect(win, width=50, height=50, pos=(100, -200.0))
    clock       = core.Clock()
    
    # show button options as colored squares in instructions
    stim_sq1.setColor(subj_colors[0],'rgb255')
    stim_sq2.setColor(subj_colors[1],'rgb255')
    stim_sq3.setColor(subj_colors[2],'rgb255')
    stim_sq4.setColor(subj_colors[3],'rgb255')
    

    # Welcome instructions
    stim_instr.setText(welcome_txt)
    stim_instr.draw()
    stim_sq1.draw() # draw button options
    stim_sq2.draw() # draw button options
    stim_sq3.draw() # draw button options
    stim_sq4.draw() # draw button options
    win.flip()
    core.wait(2)
    event.waitKeys(keyList = buttons+['q']) # MONITOR FOR POSSIBLE KEY PRESS
                
     # q quits the experiment
    if button_pressed == 'q':
        core.quit()
        
        # For quitting early
        keys = event.getKeys()
        if keys:
            # q quits the experiment
            if keys[0] == 'q':
                core.quit()
        
# Close-up   
win.close() # Close window
core.quit()














