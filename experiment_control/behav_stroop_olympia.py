"""
Single-letter behavioral stroop task
O.Colizoli 2023
"""
# data saved in ~/logfiles_main/sub-XXX

# Import necessary modules
import os, time  # for paths and data
import random
from psychopy import core, visual, event, data, sound, gui, monitors
import numpy as np
import pandas as pd
import general_parameters as gp # for letter sets, counterbalancing
#from rusocsci import buttonbox
#from IPython import embed as shell

debug_mode = False

sqw = 150 # colored square width
sqh = 150 # colored square height
fh  = 50  # fixation cross height
ww = 1000 # wrap width of instructions text

### PARAMETERS ###
# Timing in seconds
REPS    = 12   # how many times to repeat all trials in training_balancing
t_black = .2   # black letter
t_stim  = 1.5  # stimulus duration of colored letter, maximum RT
t_ITI   = [1,2.5]  # feedback (ITI)
        
## load CSV file for balancing stimuli
stroop_trials = pd.read_csv(os.path.join('stimuli','stroop_balancing.csv'))

## Buttons 
response_buttons = ["z","x","n","m"]

# Get subject number
g = gui.Dlg()
g.addField('Subject Number:')
g.show()
subject_ID = int(g.data[0])

## Create LogFile folder cwd/LogFiles
cwd = os.getcwd()
logfile_dir = os.path.join(cwd,'logfiles_main','sub-{}'.format(subject_ID)) 
if not os.path.isdir(logfile_dir):
    os.makedirs(logfile_dir)

if subject_ID:
    ### Define stimuli for this subject & save file
    subject_colors  = os.path.join(cwd,'colors','sub-{}'.format(subject_ID),'sub-{}_colors.csv'.format(subject_ID))
    if not os.path.exists(subject_colors): # only create stimuli if file does not exist for this subject
        break
    else: # load existing file
        DFS = pd.read_csv(subject_colors)
        subj_colors = [[DFS['r'][i],DFS['g'][i],DFS['b'][i]] for i in range(len(DFS))]
        subj_letters = list(DFS['letter'])
    
#    if np.mod(subject_ID,2) == 0:
#        buttons = p.buttons # 0,1
#        button_names = p.button_names
#    else:
#        buttons = [p.buttons[1],p.buttons[0]]
#        button_names = [p.button_names[1],p.button_names[0]]
    
    # Set-up button box
    #bb = buttonbox.Buttonbox()    
    
    ## output file name with time stamp prevents any overwriting of data
    timestr = time.strftime("%Y%m%d-%H%M%S") 
    output_filename = os.path.join(logfile_dir,'sub-{}_task-stroop_{}.csv'.format(subject_ID,timestr ))
    # output dataframe
    cols = ['subject','session','group','trial_num','letter','r','g','b','congruent','button_pressed','correct','RT','ITI']
    DF = pd.DataFrame(columns=cols)
    
    # Set-up window:
    mon = monitors.Monitor('myMac15', width=p.screen_width, distance=p.screen_dist)
    mon.setSizePix((p.scnWidth, p.scnHeight))
    win = visual.Window((p.scnWidth, p.scnHeight),color=p.grey,colorSpace='rgb255',monitor=mon,fullscr=not debug_mode,units='pix',allowStencil=True,autoLog=False)
    win.setMouseVisible(False)
    
    # Set-up stimuli and timing
    welcome_txt = 'Single-letter Stroop Task! \
    \nYou will see letters appear one at a time on the screen.\
    \nInstructions: Indicate the COLOR of the letter as fast and as accurately as possible.\
    \n\n\n<Press any button to continue>'
    
    pauze = 'Take a short break now.\
    \n\n\n<Press any button to CONTINUE the task>'

    stim_instr  = visual.TextStim(win, color='black', pos=(0.0, 0.0), wrapWidth=ww)  # can't really center, and TextBox doesn't work, stupid!
    # odd_sq      = visual.Rect(win, width=sqw,height=sqh,autoLog=None,pos=(300, 25))
    stim_fix    = visual.TextStim(win, text='+', color='black', pos=(0.0, 0.0), height=fh)
    # stim_sq     = visual.Rect(win, width=sqw, height=sqh, pos=(0.0, 0.0))
    feed_miss   = visual.TextStim(win, text='Too slow!', color='blue', pos=(0.0, 50), wrapWidth=8)  # can't really center, and TextBox doesn't work, stupid!
    clock       = core.Clock()
    
    # Set conditions and stimulus list    
    ALL_TRIALS  = pd.concat([stroop_trials]*REPS, ignore_index=True)    
    ALL_TRIALS = ALL_TRIALS.sample(frac=1).reset_index(drop=True) ## SHUFFLE ORDER OF TRIALS
    if debug_mode:
        ALL_TRIALS = ALL_TRIALS.iloc[0:5,:] # just get first N trials

    # Welcome instructions
    stim_instr.setText(welcome_txt)
    stim_instr.draw()
    win.flip()
    core.wait(0.25)
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
                win.flip()
                core.wait(0.5)
                #b = bb.waitButtons()
                # Wait a few seconds before first trial to stabilize gaze
                stim_fix.setColor('black')
                stim_fix.draw()
                win.flip()
                core.wait(3)
        except:
            pass
        
        stim_fix.setColor('black')
        # current trial values
        this_color_idx  = int(ALL_TRIALS['color'][t]) 
        this_letter_idx = int(ALL_TRIALS['letter'][t])
        this_freq = ALL_TRIALS['frequency'][t]
        this_actual = ALL_TRIALS['actual_frequency'][t]
        this_answer = ALL_TRIALS['oddball'][t]
        
        # Check if oddball LETTER or NUMBER?
        this_letter = []
        if this_letter_idx == 7: # oddball number
            this_letter = oddball_numbers[random.randint(0, 8)]  # choose random 1-9
        elif this_letter_idx == 100: # oddball color + letter
            this_letter = subj_letters[random.randint(0, 3)] # choose random letter 1-4
        else:
            this_letter = subj_letters[this_letter_idx-1]
            
        # Check if oddball COLOR?
        this_color = []
        if this_color_idx == 7: # odd ball color!!
            this_color = subj_oddcolor
        elif this_color_idx == 100: # random normal color, always paired with a number
            this_color_idx = random.randint(1, 4)
            this_color = subj_colors[this_color_idx-1] # index 0-3, not 1-4
        else:
            this_color = subj_colors[this_color_idx-1]
            
        # Target stimuli current trial        
        print('########## Trial {} #########'.format(t+1))
        print('Color of current trial: {}'.format(this_color))
        print('Letter of current trial: {}'.format(this_letter))
        print('Correct answer is: {}'.format(buttons[this_answer]))
        
        # STIMULUS (letter) & wait for response
        stim_letter.setSound(os.path.join('stimuli','{}.wav'.format(this_letter)))
        # stim_sq.setImage(os.path.join('stimuli','{}.png'.format(this_color)))
        stim_sq.lineColorSpace = 'rgb255'
        stim_sq.lineColor = p.grey # grey
        stim_sq.fillColorSpace = 'rgb255'
        stim_sq.fillColor = this_color

        respond1 = [] # during stim
        respond2 = [] # after stim
        bb.clearEvents() # clear button box events

        # STIMULUS image + audio letter
        stim_letter.play() 
        core.wait(0.5) # syncs sound onset with visual onset
        clock.reset() # for latency measurements
        stim_sq.draw()    
        win.flip()

        # check if response made during stimulus
        respond1 = bb.waitButtons(maxWait=t_stim, buttonList=buttons, timeStamped=clock)
        
        # 2nd RESPONSE INTERVAL fixation cross
        if not respond1:
            # did not respond during stimulus, see if they will respond during fixation:
            # show fixation cross
            stim_fix.setColor('black')
            stim_fix.draw() 
            win.flip()     
            respond2 = bb.waitButtons(maxWait=t_RT, buttonList=buttons, timeStamped=clock)
        else:
            # responded during stimulus, keep stimulus on screen
            response, latency = respond1[0]
            RT = round(latency,8)
            core.wait(t_stim-latency) # wait for the rest of the stimulus period
          
        ## get RT and ACCURACY
        if respond2: # responded after stimulus
            response, latency = respond2[0]
            RT = round(latency,8)
            core.wait(t_RT-RT) # wait for the rest of the ITI period
        elif respond1: # responded during stimulus
            stim_fix.setColor('black')
            stim_fix.draw() 
            win.flip()   
            core.wait(t_RT) # wait for the entire ITI period
        else: # missed
            stim_fix.setColor('black')
            stim_fix.draw() 
            win.flip()   
            response, RT = ('missing', np.nan)

        ## check whether correct response?
        correct = buttons[this_answer] == response
                
        ## FEEDBACK
        if (response == 'missing') or (not correct):
            if response == 'missing':
                feed_miss.draw()
                stim_fix.setColor('blue')
                stim_fix.draw() 
            else:
                feed_error.draw()
                stim_fix.setColor('red')
                stim_fix.draw() 
        else:
            feed_correct.draw()
            stim_fix.setColor('green')
            stim_fix.draw() 
        # ITI - fixation with feedback
        win.flip()
        core.wait(t_ITI)

        # output data frame on each trial 
        DF.loc[t] = [
            int(subject_ID),    # subject
            int(t),             # trial_num
            this_letter,        # letter
            int(this_freq),     # frequency
            int(this_actual),   # actual frequency
            #this_color,         # image
            int(this_color[0]), # r
            int(this_color[1]), # g
            int(this_color[2]), # b
            int(this_answer),   # correct answer, 1=oddball,0=normal
            response,           # button pressed
            correct,            # correct response or erro
            RT,                 # RT
            t_ITI,              # ITI
           ]
        DF.to_csv(output_filename)
        
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














