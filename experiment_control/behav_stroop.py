"""
Single-letter behavioral stroop task
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
# from IPython import embed as shell

debug_mode = False

cwd = os.getcwd() # current working directory

stroop_letter_conditions = [
    ['e', 'n', 'o', 'r'], # counter balanced letters
    ['a', 'i', 's', 't' ]
]

## Keyboard Buttons 
buttons = ["z","x","c","v"] # colors 1,2,3,4 in order

### PARAMETERS ###
# Timing in seconds
REPS    = 10   # how many times to repeat all trials in trials_stroop
t_black = .2   # black letter
t_stim  = 1.5  # stimulus duration of colored letter, maximum RT
t_ITI   = [1,1.5]  # feedback (ITI)
        
## load CSV file for balancing stimuli
stroop_trials = pd.read_csv(os.path.join(cwd, 'stimuli','trials_stroop.csv'),sep='\t')

# Get subject number
g = gui.Dlg()
g.addField('Subject Number:')
g.addField('Session:')
g.show()
subject_ID = int(g.data[0])
session = int(g.data[1])
group = [3 if subject_ID < 400 else 4]

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
    
    # Set-up button box
    # bb = buttonbox.Buttonbox()    
    
    ## Create output directory cwd/sub-xxx/ses-x/logfile
    logfile_dir = os.path.join(cwd,'source','sub-{}'.format(subject_ID),'ses-{}'.format(session)) 
    if not os.path.isdir(logfile_dir):
        os.makedirs(logfile_dir)
    ## output file name with time stamp prevents any overwriting of data
    timestr = time.strftime("%Y%m%d-%H%M%S") 
    output_filename = os.path.join(logfile_dir,'sub-{}_ses-{}_task-stroop_{}.csv'.format(subject_ID,session,timestr ))
    # output dataframe
    cols = ['subject','session','group','trial_num','letter','r','g','b','congruent','correct_button','button_pressed','correct','RT','ITI','black_letter_duration']
    DF = pd.DataFrame(columns=cols)
    
    # Set-up window:
    mon = monitors.Monitor('myMac15', width=p.screen_width, distance=p.screen_dist)
    mon.setSizePix((p.scnWidth, p.scnHeight))
    win = visual.Window((p.scnWidth, p.scnHeight),color=p.white,colorSpace='rgb255',monitor=mon,fullscr=not debug_mode,units='pix',allowStencil=True,autoLog=False)
    win.setMouseVisible(False)
    
    # Set-up stimuli and timing
    welcome_txt = 'Single-Letter Stroop Task! \
    \n\nYou will see letters appear one at a time on the screen.\
    \nInstructions: Indicate the COLOR of the letter as fast AND as accurately as possible.\
    \nThe response options are shown in the colored squares below.\
    \nPlease use only one hand to respond - your dominant hand.\
    \n\n\n<Press any button to continue>'
    
    pauze = 'Take a short break now.\
    \n\n\n<Press any button to CONTINUE the task>'

    stim_instr  = visual.TextStim(win, color='black', pos=(0.0, 0.0), wrapWidth=p.ww)  # can't really center, and TextBox doesn't work, stupid!
    stim_fix    = visual.TextStim(win, text='+', color='black', pos=(0.0, 0.0), height=p.fh)
    stim_letter = visual.TextStim(win, color ='black', pos=(0.0, p.sp), height=p.lh2, font=p.font_trained) 
    feed_miss   = visual.TextStim(win, text='Too slow!', color='blue', pos=(0.0, 50.0))  # can't really center, and TextBox doesn't work, stupid!
    stim_sq1    = visual.Rect(win, width=50, height=50, pos=(-100, -200))
    stim_sq2    = visual.Rect(win, width=50, height=50, pos=(-50, -200))
    stim_sq3    = visual.Rect(win, width=50, height=50, pos=(0, -200))
    stim_sq4    = visual.Rect(win, width=50, height=50, pos=(50, -200))
    clock       = core.Clock()
    # show button options as colored squares
    stim_sq1.setColor(subj_colors[0],'rgb255')
    stim_sq2.setColor(subj_colors[1],'rgb255')
    stim_sq3.setColor(subj_colors[2],'rgb255')
    stim_sq4.setColor(subj_colors[3],'rgb255')
    
    # Set conditions and stimulus list    
    ALL_TRIALS  = pd.concat([stroop_trials]*REPS, ignore_index=True)    
    ALL_TRIALS = ALL_TRIALS.sample(frac=1).reset_index(drop=True) ## SHUFFLE ORDER OF TRIALS
    # if debug_mode:
    #     ALL_TRIALS = ALL_TRIALS.iloc[0:5,:] # just get first N trials
        
    # Welcome instructions
    stim_instr.setText(welcome_txt)
    stim_instr.draw()
    stim_sq1.draw() # draw button options
    stim_sq2.draw() # draw button options
    stim_sq3.draw() # draw button options
    stim_sq4.draw() # draw button options
    win.flip()
    core.wait(1)
    event.waitKeys(keyList = buttons+['q']) # MONITOR FOR POSSIBLE KEY PRESS
    
    #b = bb.waitButtons()
    
    # Wait a few seconds before first trial
    stim_fix.draw()
    win.flip()
    core.wait(3)
    
    #### TRIAL LOOP ###
    for t in range(len(ALL_TRIALS)):
        
        try:
            if (t==np.floor(len(ALL_TRIALS)/3)) or (t==np.floor((len(ALL_TRIALS)/3)*2)):
                # take a break!
                stim_instr.setText(pauze)
                stim_instr.draw()
                stim_sq1.draw() # draw button options
                stim_sq2.draw() # draw button options
                stim_sq3.draw() # draw button options
                stim_sq4.draw() # draw button options
                win.flip()
                core.wait(0.5)
                event.waitKeys(keyList = buttons+['q']) # MONITOR FOR POSSIBLE KEY PRESS
                #b = bb.waitButtons()
                # Wait a few seconds before first trial to stabilize gaze
                stim_fix.setColor('black')
                stim_fix.draw()
                win.flip()
                core.wait(3)
        except:
            pass
                
        # current trial letter + color values
        this_color_idx  = int(ALL_TRIALS['color'][t]) - 1
        this_letter_idx = int(ALL_TRIALS['letter'][t]) - 1
        
        # get current letter
        this_letter = subj_letters[this_letter_idx]
        
        # get current color
        this_color = subj_colors[this_color_idx]
        
        # get correct button
        correct_button = buttons[this_color_idx]
            
        # Target stimuli current trial        
        print('########## Trial {} #########'.format(t+1))
        print('Color of current trial: {}'.format(this_color))
        print('Letter of current trial: {}'.format(this_letter))
        
        # BLACK LETTER STIMULUS 200 ms fixed duration
        stim_letter.setText(this_letter)
        stim_letter.setColor('black')
        stim_letter.draw()
        win.flip()
        stim_onset = clock.getTime() # letter stimulus onset locked to black letter
        core.wait(t_black)
        
        # COLORED LETTER STIMULUS wait for response
        stim_letter.setText(this_letter)
        stim_letter.setColor(this_color,'rgb255')
        stim_letter.draw()
        win.flip()
        color_onset = clock.getTime()
        # response = event.waitKeys(keyList=buttons+['q'],timeStamped=clock, maxWait=t_stim) # Maximum duration
        response = event.waitKeys(keyList=buttons+['q'],timeStamped=clock) # No maximum duration
        
        # if too slow/missed, warn the participant!
        if not response:
            button_pressed = np.nan   # didn't press
            RT = np.nan
            # warning!
            feed_miss.draw()
            stim_fix.setColor('blue')
            stim_fix.draw()
            win.flip()
            core.wait(1.5)
        else: # button pressed + RT
            button_pressed = response[0][0]
            RT = response[0][1] - stim_onset # clock time minus stim onset of black letter
            
        # q quits the experiment
        if button_pressed == 'q':
            core.quit()
        
        # Inter-trial interval, jittered
        ITI = np.round(random.uniform(t_ITI[0],t_ITI[1]),2) # draw from uniform distribution
        stim_fix.setColor('black')
        stim_fix.draw()
        win.flip()
        core.wait(ITI)
        
        # output data frame on each trial, duration (seconds)  
        # cols = ['subject','session','group','trial_num','letter','r','g','b','congruent','correct_button','button_pressed','correct','RT','ITI','black_letter_duration']
        DF.loc[t] = [
            int(subject_ID),                    # subject
            session,                            # session
            group[0],                           # group
            t,                                  # trial_number
            this_letter,                        # letter
            int(this_color[0]),                 # r
            int(this_color[1]),                 # g
            int(this_color[2]),                 # b
            int(ALL_TRIALS['congruent'][t]),    # congruent boolean
            correct_button,                     # correct button
            button_pressed,                     # button pressed by participant
            correct_button == button_pressed,   # correct response if correct_button == button_pressed
            round(RT,8),                        # RT
            ITI,                                # ITI
            color_onset - stim_onset            # black letter duration
                    ]     
                                    
        DF.to_csv(output_filename,sep='\t')
        
        # For quitting early
        keys = event.getKeys()
        if keys:
            # q quits the experiment
            if keys[0] == 'q':
                core.quit()
           
    # End screen for participants
    stim_instr.setText('Well done!')
    stim_instr.draw()
    win.flip()
    core.wait(2)
        
# Close-up   
win.close() # Close window
core.quit()














