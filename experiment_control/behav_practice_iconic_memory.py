#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
ICONIC MEMORY TASK
O.Colizoli & M.Blasco Oliver, 2019
"""
# no data saved for practice trials

# Import necessary modules
from psychopy import core, visual, event, gui, monitors, data, misc
import random
import numpy as np
import os, time # for paths and data
import pandas as pd
from IPython import embed as shell
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

## see general_parameters.py
cue             = gp.cue
pos_pol         = gp.pos_pol

#Sizes
letter_size     = gp.letter_size
eccen           = gp.eccen 
cue_size        = gp.cue_size 
text_size       = gp.text_size
fixcross_size   = gp.fixcross_size
y_shift_up      = gp.y_shift_up

# INTEGRATE DISTANCE WITH THE MONITOR SOMWHERE: 73 cm according to paper
#Other colors:
cue_target_color    = gp.cue_target_color
cue_holder_color    = gp.cue_holder_color
background_color    = gp.background_color
instructions_color  = gp.instructions_color
letters_color       = gp.letters_color

# Timing
read_time             = gp.read_time
fix_time              = gp.fix_time
letters_time          = gp.letters_time
cue_time              = gp.cue_time
max_response_duration = gp.max_response_duration
endtrial_time         = gp.endtrial_time

# Response buttons values
square_fillcolor = gp.square_fillcolor
square_vertices  = gp.square_vertices
square_positions = gp.square_positions

######################################################
###################### TEXTS ########################
######################################################
welcome_txt = "PRACTICE Sensory Memory Test\
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
\nBut first you will do a practice round. \
\nIf you have questions, please let us know. Otherwise, feel free to start. Good luck!\
\n\n[Press the SPACEBAR to begin the PRACTICE]"

end_txt = "Well done! You finished the practice round. Please, call the experimenter."

######################################################
#################### EXPERIMENT #####################
######################################################

if subject_ID:
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
    # One type, numbers practice
    block_letter_list = gp.letter_sets_IM_PRACTICE
    
    # Set conditions and stimulus list    
    BTRIALS  = pd.concat([trials_IM]*REPS, ignore_index=True)
    BTRIALS = BTRIALS.sample(frac=1).reset_index(drop=True) ## SHUFFLE ORDER OF TRIALS
    BTRIALS = BTRIALS.iloc[0:10,:] # just get first X trials
            
    #### TRIALS WITHIN BLOCK LOOP ###
    for t in range(len(BTRIALS)):
        
        # Set invisible mouse
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
        if t == 0:
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

        # ITI between letters disappearing and cues
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
        # win.flip()
        # core.wait(cue_time) #show the cues for 2 seconds
        
        # then show the 8 option to select with mouse
        # myMouse.setPos([0,-8.75])
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
                break
                                   
        # Blank screen
        fixcross.setAutoDraw(False)
        win.flip()
        core.wait(endtrial_time)
        
    # End screen for participants
    instr.setText(end_txt)
    instr.draw()
    win.flip()
    event.waitKeys()


# Close-up  
win.close()
core.quit()














