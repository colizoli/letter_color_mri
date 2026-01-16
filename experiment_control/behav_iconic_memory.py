#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
ICONIC MEMORY TASK
O.Colizoli & M.Blasco Oliver, 2019
"""
# data saved in ~/source/sub-XXX/ses-X

# Import necessary modules
from psychopy import core, visual, event, gui, monitors, data, misc
import random
import numpy as np
import os, time # for paths and data
import pandas as pd
#from IPython import embed as shell
import general_parameters as gp # for letter sets, counterbalancing, iconic memory parameters

######################################################
##################### PARAMETERS #####################
######################################################
debug_mode = False

## load CSV file for balancing stimuli
trials_IM = pd.read_csv(os.path.join('stimuli','trials_iconic_memory.csv'))
REPS = 1 # number of times to repeat trials_IM

#Get subject number
g = gui.Dlg()
g.addField('Subject Number:')
g.addField('Session:')
g.show()
subject_ID = int(g.data[0])
session = int(g.data[1])
# in case GUI doesn't work
# subject_ID = 1
# session = 1

break_trials    = [48,48*2,48*3,192+48,192+48*2,192+48*3]
counterbalance  = [0,0,1,1]*1000 # counter balance order letter sets, also alternate with trained/untrained sets
# Stimuli

# ICONIC MEMORY TASK PARAMETERS
cue             = '*' # red asterisk, should it be uncolored?
pos_pol         = [0,45,90,135,180,225,270,315] # positions in degrees

#Sizes
letter_size     = 1.75# 1.75 degrees in code correspond to 1.15-1.17cm on the screen which are real 1.2 degrees (like paper)
eccen           = 5.2 # 5.2 degrees in code correspond to 5 cm on the screen which are real 5.2 degrees (like paper)
cue_size        = 0.6 # 0.6 degrees in code correspond to 0.18-0.2 cm on the screen which are real 0.2 degrees (like paper 0.1 x2 cause too small)
text_size       = 28 # in pix
fixcross_size   = 2.1 # 2.1 degrees in code correspond to 0.97 cm on the screen which are real 1 degrees
y_shift_up      = 2

# INTEGRATE DISTANCE WITH THE MONITOR SOMWHERE: 73 cm according to paper
#Other colors:
cue_target_color    = gp.black  #'black' # target cue in red, holder cue in black
cue_holder_color    = [200,200,200] #gp.grey #'grey'
background_color    = gp.white #'white' # code??
instructions_color  = gp.black #'grey'
letters_color       = gp.black # all black

# Timing
read_time             = 3
fix_time              = 1 #1000ms
letters_time          = 0.1# 0.1#0.500#0.100 #100ms, presentation of the letters array
cue_time              = 1.5 #Cues remain 
max_response_duration = 4 #they click on the letter they think it was cued
endtrial_time         = 1.5 #0.05 #50ms, fix cross disappears here. This is the ITI that was 0.8 and from piloting we want to make slower. Try 1.5 s?

square_fillcolor    = cue_holder_color
square_vertices     = ((-30, -30), (30, -30), (30, 30),(-30,30))
square_positions    = [(-150,-350), (-50,-350), (50,-350), (150,-350),
                    (-150,-450),(-50,-450),(50,-450),(150,-450)]

######################################################
###################### TEXTS ########################
######################################################
welcome_txt = "Sensory Memory Test\
\nThis task is about remembering which letter was presented at which location of the screen.\
\n\nAt the beginning of each trial, you see a cross [+] in the center of the screen, always maintain fixation there.\
\nThen, 8 letters arranged on a circle are presented for a very brief time. Please, remain fixating your eyes at the [+] while this happens!\
\nOnce the letters disappear, one of the 8 letter locations is indicated with a small black asterisk [*].\
\nThe other locations are marked with a grey asterisk (ignore them).\
\n\nYour task is to remember which letter was presented at the location cued with a black asterisk.\
\nTo give your answer, a letter panel appears in the lower part of the screen.\
\nYou can use the mouse to click on the letter you think was presented in the cued location.\
\nAfter a short delay, another trial starts.\
\n\nThere is a total of 8 blocks. Each block lasts around 3-4 minutes, and you can take breaks in between.\
\nIf you have questions, please let us know. Otherwise, feel free to start. Good luck!\
\n\n[Press the SPACEBAR to begin]"

break_txt = "Well done! You have completed {}/8 blocks. \n\
You can take a short break now. \n\n\n\n\
*** Press [SPACEBAR] to continue ***"

break_block_txt = "Well done! You have completed {}/8 blocks. \n\
You can take a short break now. \n\n\n\
Note that from now on, the letters you will see will be different.\n\n\n\
*** Press [SPACEBAR] to continue ***"

end_txt = 'Well done! You have finished the experiment.\nThank you!\n\n\n\n\
*** Press [SPACEBAR] to EXIT the experiment ***'

######################################################
#################### EXPERIMENT #####################
######################################################

if subject_ID:
    
    ## Create LogFile folder cwd/LogFiles
    cwd = os.getcwd()
    # bids format    
    logfile_dir = os.path.join(cwd,'source','sub-{}'.format(subject_ID),'ses-{}'.format(session),'behav')  # system independent
    
    if not os.path.isdir(logfile_dir):
        os.makedirs(logfile_dir)
        
    ## output file name with time stamp prevents any overwriting of data
    timestr = time.strftime("%Y%m%d-%H%M%S") 
    output_filename = os.path.join(logfile_dir,'sub-{}_ses-{}_task-iconic_behav_{}.tsv'.format(subject_ID,session,timestr))
    header = ['subject',
              'session',
              'letter_set',
              'trained',
              'trial_number',
              'block_trial_number',
              'target_position',
              'target_letter',
              'ISI',
              'button_clicked',
              'RT',
              'correct',
              'position_list',
              'letter_list'
          ]
    DF = pd.DataFrame(columns=header) #called DF in other script

    #Set-up window:    
    mon = monitors.Monitor('behavlab', width=gp.screen_width, distance=gp.screen_dist)
    mon.setSizePix((gp.scnWidth, gp.scnHeight))
    win = visual.Window((gp.scnWidth, gp.scnHeight),color=background_color,colorSpace='rgb255',monitor=mon,fullscr=not debug_mode,units='deg')
    win.update()
    myMouse = event.Mouse(visible = True, win = win)

    #Set-up stimuli and timing
    fixcross        = visual.TextStim(win, text='+', colorSpace = 'rgb255', color=[0,0,0], units = "deg", pos = [0,y_shift_up], height=fixcross_size)
    instr           = visual.TextStim(win, colorSpace = 'rgb255', color=instructions_color, units = 'pix', height=text_size,wrapWidth=gp.ww)
    stim_cue        = visual.TextStim(win, colorSpace = 'rgb255', text=cue)
    stim_letters    = visual.TextStim(win, font = gp.font_trained, colorSpace = 'rgb255', color=[0,0,0], height=letter_size)
    clock           = core.Clock() 
    
    # Response options with mouse:
    square_1 = visual.ShapeStim(win, units='pix',lineWidth=2, lineColor='black',fillColor=square_fillcolor, fillColorSpace='rgb255',
                               vertices=square_vertices, closeShape=True,
                               pos=square_positions[0])
    square_2 = visual.ShapeStim(win, units='pix',lineWidth=2, lineColor='black',fillColor=square_fillcolor, fillColorSpace='rgb255',
                               vertices=square_vertices, closeShape=True,
                               pos=square_positions[1])
    square_3 = visual.ShapeStim(win, units='pix',lineWidth=2, lineColor='black',fillColor=square_fillcolor, fillColorSpace='rgb255',
                               vertices=square_vertices, closeShape=True,
                               pos=square_positions[2])
    square_4 = visual.ShapeStim(win, units='pix',lineWidth=2, lineColor='black',fillColor=square_fillcolor, fillColorSpace='rgb255',
                               vertices=square_vertices, closeShape=True,
                               pos=square_positions[3])
    square_5 = visual.ShapeStim(win, units='pix',lineWidth=2, lineColor='black',fillColor=square_fillcolor, fillColorSpace='rgb255',
                               vertices=square_vertices, closeShape=True,
                               pos=square_positions[4])
    square_6 = visual.ShapeStim(win, units='pix',lineWidth=2, lineColor='black',fillColor=square_fillcolor, fillColorSpace='rgb255',
                               vertices=square_vertices, closeShape=True,
                               pos=square_positions[5])
    square_7 = visual.ShapeStim(win, units='pix',lineWidth=2, lineColor='black',fillColor=square_fillcolor, fillColorSpace='rgb255',
                               vertices=square_vertices, closeShape=True,
                               pos=square_positions[6])
    square_8 = visual.ShapeStim(win, units='pix',lineWidth=2, lineColor='black',fillColor=square_fillcolor, fillColorSpace='rgb255',
                               vertices=square_vertices, closeShape=True,
                               pos=square_positions[7])
    squares_stimholder = [square_1, square_2, square_3, square_4, square_5, square_6, square_7, square_8]
    
    
    letter_1 = visual.TextStim(win,units='pix',font = gp.font_trained, height = 30, color='black',pos=square_positions[0])
    letter_2 = visual.TextStim(win,units='pix',font = gp.font_trained, height = 30, color='black',pos=square_positions[1])
    letter_3 = visual.TextStim(win,units='pix',font = gp.font_trained, height = 30, color='black',pos=square_positions[2])
    letter_4 = visual.TextStim(win,units='pix',font = gp.font_trained, height = 30, color='black',pos=square_positions[3])
    letter_5 = visual.TextStim(win,units='pix',font = gp.font_trained, height = 30, color='black',pos=square_positions[4])
    letter_6 = visual.TextStim(win,units='pix',font = gp.font_trained, height = 30, color='black',pos=square_positions[5])
    letter_7 = visual.TextStim(win,units='pix',font = gp.font_trained, height = 30, color='black',pos=square_positions[6])
    letter_8 = visual.TextStim(win,units='pix',font = gp.font_trained, height = 30, color='black',pos=square_positions[7])
    
    letters_stimholder = [letter_1,letter_2,letter_3,letter_4,letter_5,letter_6,letter_7,letter_8]    
    
    # Welcome instructions
    instr.setText(welcome_txt)
    instr.draw()
    win.flip()
    event.waitKeys()
    
    #### MAIN LOOP ###
    # Two main blocks of Trained and Untrained
    trial_counter = 0
    breaks_counter = 1
    blocknumber = 1
    block_order = [int(counterbalance[subject_ID]), int(not counterbalance[subject_ID])]
    
    for letter_set in block_order:
        
        block_letter_list = gp.letter_sets_lower_8_IM[letter_set] # COUNTER BALANCED
        
        # Set conditions and stimulus list    
        BTRIALS  = pd.concat([trials_IM]*REPS, ignore_index=True)
        BTRIALS = BTRIALS.sample(frac=1).reset_index(drop=True) ## SHUFFLE ORDER OF TRIALS
        if debug_mode:
            BTRIALS = BTRIALS.iloc[0:3,:] # just get first X trials
        
        # Break between blocks
        print(blocknumber)
        if blocknumber > 1:
            instr.setText(break_block_txt.format(breaks_counter))
            instr.draw()
            win.flip()
            core.wait(0.001)
            event.waitKeys()
            print('Break: ', breaks_counter)
            breaks_counter +=1
            #print('Break: ', breaks_counter)

        #### TRIALS WITHIN BLOCK LOOP ###
        for t in range(len(BTRIALS)):   
            
            if trial_counter in break_trials:
                instr.setText(break_txt.format(breaks_counter))
                instr.draw()
                win.flip()
                event.waitKeys()
                print('Break: ', breaks_counter)
                breaks_counter +=1
                #print('Break: ', breaks_counter)
            
            myMouse.setVisible(0)
            
            # Target stimuli current trial
            target_letter   = block_letter_list[int(BTRIALS['letter'][t])]
            target_pos_pol  = BTRIALS['position'][t]
            target_ISI      = BTRIALS['ISI'][t]
        
            # Generate non-target positions and letters, and shuffle orders
            pos_list = list(pos_pol) # all possible positions
            pos_list.remove(target_pos_pol)
            random.shuffle(pos_list)

            letter_list = list(block_letter_list) # all possible letters
            letter_list.remove(target_letter)
            random.shuffle(letter_list)
       
            # Fixation cross
            fixcross.setAutoDraw(True)
            win.flip()
            if trial_counter in break_trials or trial_counter == 0:
                core.wait(fix_time+1.5)
            else:
                core.wait(fix_time) #0.3 Wait 300 ms so that they fixate, then show letters
                            
            #Present 8 letters
            for l in range(len(block_letter_list)): #just to have 8
                if l == 0:
                    current_pos = misc.pol2cart(target_pos_pol,eccen,units='deg')
                    stim_letters.setPos([current_pos[0],current_pos[1]+y_shift_up])
                    #stim_letters.setPos(misc.pol2cart(target_pos_pol,eccen,units='deg')) # change
                    stim_letters.setText(target_letter)
                else:
                    current_pos = misc.pol2cart(pos_list[l-1],eccen,units='deg')
                    stim_letters.setPos([current_pos[0],current_pos[1]+y_shift_up]) # change
                    stim_letters.setText(letter_list[l-1][0]) #letters from pair
                stim_letters.draw() 
            win.flip()
            
            #Present letters for certain ms and then disappear
            core.wait(letters_time) #now longer to see
            win.flip()
            
            # ISI between letters disappearing and cues
            core.wait(target_ISI)
        
            #Present cues (big/red or small/black depending on target location)
            for c in range(len(block_letter_list)):
                if c == 0: #Target position-letter
                    stim_cue.setHeight(cue_size)
                    current_pos = misc.pol2cart(target_pos_pol,eccen,units='deg')
                    stim_cue.setPos([current_pos[0],current_pos[1]+y_shift_up])
                    stim_cue.setColor(cue_target_color)
                else: #Other positions-letters
                    stim_cue.setHeight(cue_size)
                    current_pos = misc.pol2cart(pos_list[c-1],eccen,units='deg')
                    stim_cue.setPos([current_pos[0],current_pos[1]+y_shift_up]) # change
                    stim_cue.setColor(cue_holder_color)
                
                stim_cue.draw()
            
            # then show the 8 option to select with mouse
            myMouse.setPos([0,-11.8]) #put the mouse in the center of the letter oannel (in deg units)
            myMouse.setVisible(1)
            for s in range(len(squares_stimholder)):
                squares_stimholder[s].draw()
                letters_stimholder[s].setText(block_letter_list[s])
                letters_stimholder[s].draw()
          
            win.flip()
            
            # get mouse response:
            onset_time = core.getTime()
            RT = None
            response = None
            while RT is None and (core.getTime() - onset_time) <= max_response_duration:  # for duration time we sample a response
            
                mouse_click = myMouse.getPressed()
                mouse_click = mouse_click[0]
                
                # For quitting earlier
                keys = event.getKeys()
                if keys:
                    # q quits the experiment
                    if keys[0] == 'q':
                        win.close()
                        core.quit()
                    
                #Record response and show it
                if mouse_click:
                    for resp in range(len(letters_stimholder)):
                        if myMouse.isPressedIn(squares_stimholder[resp]): #, buttons=[0]):
                            print('Clicked letter: {}'.format(block_letter_list[resp]))
                            response = block_letter_list[resp]
                            mouse_pos = myMouse.getPos() #in degrees (window units)
                            RT = core.getTime() - onset_time  # Reaction time
                            print(RT)                    
                    #break
                                      
            # Blank screen
            fixcross.setAutoDraw(False)
            win.flip()
            core.wait(endtrial_time)
            
            #After 3 seconds, no response?
            if response == None:
                response = np.nan
        
            # Accuracy
            if response == target_letter:
                correct = 1
            else:
                correct = 0
            print('Response letter: ', response)
            print('Acurracy: ', correct)            

            # Save trial data
            DF.loc[trial_counter] = [
                subject_ID,      # subject number
                session,         # Session numb er 1 or 2 (pre-post)
                letter_set,      # set of trained letters 0 or 1
                int(letter_set == np.mod(subject_ID,2)), # letter set trained or untrained?
                trial_counter,   # trial number all
                t,               # trial number within the block
                target_pos_pol,  # Target position in xy
                target_letter,   # Target letter
                target_ISI,      # Current ISI
                response,        # whether they responded d or j
                RT,              # RT of mouse answer
                correct,         # Accuracy (1 correct, 0 incorrect)
                pos_list,
                letter_list
            ]          
            DF.to_csv(output_filename,sep='\t')
            trial_counter += 1
        
        blocknumber += 1

    # End screen for participants
    instr.setText(end_txt)
    instr.draw()
    win.flip()
    event.waitKeys()

# Close-up  
win.close()
core.quit()














