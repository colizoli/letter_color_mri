#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
CHOOSING LETTER-COLOR ASSOCIATIONS TASK
O.Colizoli & M.Blasco Oliver, 2019
Outputs a TSV file for sub-101 and sub-201 (Group2 gets Group1's preferences)
"""

### Import Libraries ###
import os, time # for paths and data
from psychopy import core, visual, event, gui, monitors
import random
import numpy as np
import pandas as pd
from IPython import embed as shell
import general_parameters as gp # letter conditions, counterbalancing

debug_mode = True

#########################################
################# TEXTS #################
#########################################
welcome_txt = "Imagine that each letter of the alphabet would be associated with a color.\
\n\nYour task is to pair each letter (left hand side) to the colors presented (right hand side).\
\nThere are NO right or wrong answers, just go with your 'gut' feeling, your preferences or (dis)likes.\
\nYou may change your answers before your final decision.\
\nHowever, you may not choose the same color for more than one letter (each letter should have a unique color).\
\n\nTo make your choice, first press the letter via the keyboard, followed by ENTER.\
\nThereafter, press the number on the keyboard corresponding to the desired color, followed by ENTER.\
\nPress ENTER again, to go to the next letter.\
\n\n[Press the SPACEBAR to begin]"

######################################################
##################### FUNCTIONS #####################
######################################################

def show_input(win,response_input,type_input,options):
    keys = []
    user_txt = None
    all_options= options + ["return","backspace","space"]
    while True:
        tmp_key = event.getKeys(all_options) # keyList = [a,b,c]
        if len(tmp_key)>0:
            if tmp_key[0] == "return":
                if keys == []: # pressed forst by accident
                    continue
                else:
                    break
            elif tmp_key[0] == "backspace": # fix backspace
                try:
                    keys.pop()
                    user_txt = ''.join(keys)
                    response_input.setText(user_txt)
                    response_input.draw()
                    win.flip()
                except:
                    continue
                user_txt = ''.join(keys)
            elif tmp_key[0] == "space":
                keys.append('')
            else:
                keys.append(tmp_key[0])
                user_txt = ''.join(keys)
                response_input.setText(user_txt)
                response_input.draw()
                win.flip()
    return user_txt

######################################################
##################### PARAMETERS #####################
######################################################

#Get subject number
g = gui.Dlg()
g.addField('Subject Number:')
g.show()
subject_ID = int(g.data[0])
session = 1
# in case GUI doesn't work
# subject_ID = 1

if subject_ID:
    # trials
    trials = 8

    # Stimuli
    ## LETTER SET DEPENDENT ON SUBJECT NUMBER
    letters = gp.letter_sets_lower[np.mod(subject_ID,2)] # lower case, TRAINED set!!
    colors = gp.trained_colors

    colorcodes  = ['1.','2.','3.','4.','5.','6.','7.','8.','9.','10.','11.','12.','13.']
    colorcode_input = [1,2,3,4,5,6,7,8,9,10,11,12,13]
    colorcode_input_str = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13']

    random.shuffle(colors) # shuffle order of colors, not letters

    # Build dictionary of colors"
    color_dict = {k: 0 for k in colorcode_input}
    for k in range(len(color_dict)):
        color_dict[k+1] = colors[k]

    # Text sizes
    cc_size = 24
    l_size = 40
    cp_size = 22 # radius, so size is 32
    choice_size = 46
    choice_col_size = 16

    ### Positions X on screen
    # Letters
    let_x1      = -500           #-500
    let_x2      = let_x1+100     #-400
    # Colors number (coln) and color patches (colp)
    coln_x1     = -100           #-100
    colp_x1     = coln_x1+40     #-60
    coln_x2     = coln_x1+150    #50
    colp_x2     = coln_x2+40     #90
    # Choices (letcol) and number color (letcoln)
    letcol_x1   = 350            #350
    letcoln_x1  = letcol_x1+25   #375
    letcol_x2   = letcol_x1+70   #420
    letcoln_x2  = letcol_x2+25   #445
    letcol_xnone = np.mean([letcol_x1,letcol_x2])

    ### Positions Y on screen
    # Letters & Colors options
    let_col_y1 = np.linspace(220,-140,7).astype(int)
    distance_y = np.abs(let_col_y1[0]-let_col_y1[1])
    start_let_col_y2 = let_col_y1[0]-(distance_y/2)
    end_let_col_y2 = let_col_y1[-1]+(distance_y/2)
    let_col_y2 = np.linspace(start_let_col_y2,end_let_col_y2,6).astype(int)
    print('Position in y letters and colors:')
    print(let_col_y1)
    print(let_col_y2)

    # Choices (letcol) and number color (letcoln)
    choice_y1 = np.linspace(180,-120,7).astype(int)
    dist_y_choice = np.abs(choice_y1[0]-choice_y1[1])
    start_choice_y2 = choice_y1[0]-(dist_y_choice/2)
    end_choice_y2 = choice_y1[-1]+(dist_y_choice/2)
    choice_y2 = np.linspace(start_choice_y2,end_choice_y2,6).astype(int)
    choice_ynone = choice_y1[-1]-(dist_y_choice/2)
    print('Position in y choices:')
    print(choice_y1)
    print(choice_y2)
    print(choice_ynone)


    ### SIGN POSITIONS
    sign_y = let_col_y1[0]+60
    ask_y = let_col_y1[-1]-60
    input_y = ask_y-40
    sign_x_let = np.mean([let_x1,let_x2])
    sign_x_col = np.mean([coln_x1,colp_x2])
    sign_x_choice = np.mean([letcol_x1,letcol_x2])

    # Choosen pairs dictionary - to be filled during the task
    pairs_dict = {k: '' for k in letters}

    ## Create LogFile folder cwd/LogFiles
    cwd = os.getcwd()
    ## output file name with time stamp prevents any overwriting of data
    timestr = time.strftime("%Y%m%d-%H%M%S") 
    
    # G1
    if subject_ID < 200: # make G1 & G2's folders
        subject_yoked = subject_ID+100
        logfile_dir = os.path.join(cwd,'colors','sub-{}'.format(subject_ID))  # system independent
        logfile_dir_g2 = os.path.join(cwd,'colors','sub-{}'.format(subject_ID+100))  # system independent
        
        if not os.path.isdir(logfile_dir):
            os.makedirs(logfile_dir)
        if not os.path.isdir(logfile_dir_g2): 
            os.makedirs(logfile_dir_g2)
            
        output_filename = os.path.join(logfile_dir,'sub-{}_colors_{}.tsv'.format(subject_ID,timestr ))
        output_filename_g2 = os.path.join(logfile_dir_g2,'sub-{}_colors_{}.tsv'.format(subject_ID+100,timestr )) # G2 gets G1's colors
    # G2
    else: # save elsewhere in 'prefs' folder
        subject_yoked = subject_ID-100
        pref_dir = os.path.join(cwd,'prefs','sub-{}'.format(subject_ID))  # system independent
        if not os.path.isdir(pref_dir):
            os.makedirs(pref_dir)
        output_filename = os.path.join(pref_dir,'sub-{}_prefs_{}.tsv'.format(subject_ID,timestr ))

    header = ['subject','subject_yoked','letter','colorcode','r','g','b']
    data_pairs = pd.DataFrame(columns=header) #called DF in other script
   
    # Set-up window:
    mon = monitors.Monitor('myMac15', width=gp.screen_width, distance=gp.screen_dist)
    mon.setSizePix((gp.scnWidth, gp.scnHeight))
    win = visual.Window((gp.scnWidth, gp.scnHeight),color=gp.white,colorSpace='rgb255',monitor=mon,fullscr=not debug_mode,units='pix',allowStencil=True,autoLog=False)
    win.setMouseVisible(True)

    # Set-up stimuli and timing
    instr         = visual.TextStim(win, color='black', pos=(0.0, 0.0), wrapWidth=gp.ww ) 

    letter_sign     = visual.TextStim(win, text = 'LETTERS:', color='black', pos=(sign_x_let,sign_y), height = 24)
    color_sign      = visual.TextStim(win, text = 'COLORS:', color='black', pos=(sign_x_col,sign_y), height = 24)
    choices_sign    = visual.TextStim(win, text = 'YOUR CHOICES:', color='black', pos=(sign_x_choice,sign_y), height = 24)
    
    letter_ask      = visual.TextStim(win, text = 'Letter?', color='black', pos=(sign_x_let,ask_y),height=20)
    color_ask       = visual.TextStim(win, text = 'Color (number)?', color='black', pos=(sign_x_col,ask_y),height=20)
    letter_input    = visual.TextStim(win, color= 'black', pos=(sign_x_let,input_y),height=30)
    color_input     = visual.TextStim(win, color= 'black', pos=(sign_x_col,input_y),height=30)

    final_choice_ask = visual.TextStim(win, text = 'Continue choosing [ENTER]\nor Final choice [y]', color='black', pos=(sign_x_choice,ask_y), height = 20)
    keep_choosing    = visual.TextStim(win, text = 'Letters with no color or repeated colors. \nPlease, continue choosing. ', color='black', pos=(sign_x_choice,ask_y-40), height = 20)
    
    # Letters to be trained
    letter_01 = visual.TextStim(win, text = letters[0], color='black', font = gp.font_trained, pos=(let_x1,let_col_y1[0]), height = l_size)
    letter_02 = visual.TextStim(win, text = letters[1], color='black', font = gp.font_trained, pos=(let_x1,let_col_y1[1]), height = l_size)
    letter_03 = visual.TextStim(win, text = letters[2], color='black', font = gp.font_trained, pos=(let_x1,let_col_y1[2]), height = l_size)
    letter_04 = visual.TextStim(win, text = letters[3], color='black', font = gp.font_trained, pos=(let_x1,let_col_y1[3]), height = l_size)
    letter_05 = visual.TextStim(win, text = letters[4], color='black', font = gp.font_trained, pos=(let_x1,let_col_y1[4]), height = l_size)
    letter_06 = visual.TextStim(win, text = letters[5], color='black', font = gp.font_trained, pos=(let_x1,let_col_y1[5]), height = l_size)
    letter_07 = visual.TextStim(win, text = letters[6], color='black', font = gp.font_trained, pos=(let_x1,let_col_y1[6]), height = l_size)
    letter_08 = visual.TextStim(win, text = letters[7], color='black', font = gp.font_trained, pos=(let_x2,let_col_y2[0]), height = l_size)
    letter_09 = visual.TextStim(win, text = letters[8], color='black', font = gp.font_trained, pos=(let_x2,let_col_y2[1]), height = l_size)
    letter_10 = visual.TextStim(win, text = letters[9], color='black', font = gp.font_trained, pos=(let_x2,let_col_y2[2]), height = l_size)
    letter_11 = visual.TextStim(win, text = letters[10], color='black', font = gp.font_trained, pos=(let_x2,let_col_y2[3]), height = l_size)
    letter_12 = visual.TextStim(win, text = letters[11], color='black', font = gp.font_trained, pos=(let_x2,let_col_y2[4]), height = l_size)
    letter_13 = visual.TextStim(win, text = letters[12], color='black', font = gp.font_trained, pos=(let_x2,let_col_y2[5]), height = l_size)
    
    letter_holder = [letter_01,letter_02,letter_03,letter_04,letter_05,letter_06,letter_07,
                     letter_08,letter_09,letter_10,letter_11,letter_12,letter_13]
    
    # Color codes 01-13
    colorcode_01 = visual.TextStim(win, text = colorcodes[0], color='black', pos=(coln_x1,let_col_y1[0]), height = cc_size)
    colorcode_02 = visual.TextStim(win, text = colorcodes[1], color='black', pos=(coln_x1,let_col_y1[1]), height = cc_size)
    colorcode_03 = visual.TextStim(win, text = colorcodes[2], color='black', pos=(coln_x1,let_col_y1[2]), height = cc_size)
    colorcode_04 = visual.TextStim(win, text = colorcodes[3], color='black', pos=(coln_x1,let_col_y1[3]), height = cc_size)
    colorcode_05 = visual.TextStim(win, text = colorcodes[4], color='black', pos=(coln_x1,let_col_y1[4]), height = cc_size)
    colorcode_06 = visual.TextStim(win, text = colorcodes[5], color='black', pos=(coln_x1,let_col_y1[5]), height = cc_size)
    colorcode_07 = visual.TextStim(win, text = colorcodes[6], color='black', pos=(coln_x1,let_col_y1[6]), height = cc_size)
    colorcode_08 = visual.TextStim(win, text = colorcodes[7], color='black', pos=(coln_x2,let_col_y2[0]), height = cc_size)
    colorcode_09 = visual.TextStim(win, text = colorcodes[8], color='black', pos=(coln_x2,let_col_y2[1]), height = cc_size)
    colorcode_10 = visual.TextStim(win, text = colorcodes[9], color='black', pos=(coln_x2,let_col_y2[2]), height = cc_size)
    colorcode_11 = visual.TextStim(win, text = colorcodes[10], color='black', pos=(coln_x2,let_col_y2[3]), height = cc_size)
    colorcode_12 = visual.TextStim(win, text = colorcodes[11], color='black', pos=(coln_x2,let_col_y2[4]), height = cc_size)
    colorcode_13 = visual.TextStim(win, text = colorcodes[12], color='black', pos=(coln_x2,let_col_y2[5]), height = cc_size)
    
    colorcode_holder = [colorcode_01,colorcode_02,colorcode_03,colorcode_04,colorcode_05,colorcode_06,colorcode_07,
                        colorcode_08,colorcode_09,colorcode_10,colorcode_11,colorcode_12,colorcode_13]
    
    # Color circles RGB tones

    colorpatch_01 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[1],fillColorSpace = 'rgb255', pos =(colp_x1,let_col_y1[0]), opacity = 1.0)
    colorpatch_02 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[2],fillColorSpace = 'rgb255', pos =(colp_x1,let_col_y1[1]), opacity = 1.0)
    colorpatch_03 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[3],fillColorSpace = 'rgb255', pos =(colp_x1,let_col_y1[2]), opacity = 1.0)
    colorpatch_04 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[4],fillColorSpace = 'rgb255', pos =(colp_x1,let_col_y1[3]), opacity = 1.0)
    colorpatch_05 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[5],fillColorSpace = 'rgb255', pos =(colp_x1,let_col_y1[4]), opacity = 1.0)
    colorpatch_06 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[6],fillColorSpace = 'rgb255', pos =(colp_x1,let_col_y1[5]), opacity = 1.0)
    colorpatch_07 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[7],fillColorSpace = 'rgb255', pos =(colp_x1,let_col_y1[6]), opacity = 1.0)
    colorpatch_08 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[8],fillColorSpace = 'rgb255', pos =(colp_x2,let_col_y2[0]), opacity = 1.0)
    colorpatch_09 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[9],fillColorSpace = 'rgb255', pos =(colp_x2,let_col_y2[1]), opacity = 1.0)
    colorpatch_10 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[10],fillColorSpace = 'rgb255', pos =(colp_x2,let_col_y2[2]), opacity = 1.0)
    colorpatch_11 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[11],fillColorSpace = 'rgb255', pos =(colp_x2,let_col_y2[3]), opacity = 1.0)
    colorpatch_12 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[12],fillColorSpace = 'rgb255', pos =(colp_x2,let_col_y2[4]), opacity = 1.0)
    colorpatch_13 = visual.Circle(win=win,radius=cp_size,fillColor=color_dict[13],fillColorSpace = 'rgb255', pos =(colp_x2,let_col_y2[5]), opacity = 1.0)

    colorpatch_holder = [colorpatch_01,colorpatch_02,colorpatch_03,colorpatch_04,colorpatch_05,colorpatch_06,colorpatch_07,
                         colorpatch_08,colorpatch_09,colorpatch_10,colorpatch_11,colorpatch_12,colorpatch_13]


    # Choices holders: colored letter
#    pos_choice = [(200,150),(200,100),(200,50),(200,0),(200,-50),(200,-100),(200,-150),
#                  (250,125),(250,75),(250,25),(250,-25),(250,-75),(250,-125)]
    globals() ['choice_{}'.format(letters[0])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x1,choice_y1[0]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[1])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x1,choice_y1[1]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[2])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x1,choice_y1[2]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[3])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x1,choice_y1[3]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[4])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x1,choice_y1[4]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[5])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x1,choice_y1[5]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[6])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x1,choice_y1[6]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[7])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x2,choice_y2[0]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[8])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x2,choice_y2[1]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[9])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x2,choice_y2[2]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[10])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x2,choice_y2[3]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[11])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x2,choice_y2[4]),height=choice_size,colorSpace='rgb255')
    globals() ['choice_{}'.format(letters[12])] = visual.TextStim(win,font = gp.font_trained, pos=(letcol_x2,choice_y2[5]),height=choice_size,colorSpace='rgb255')
    choice_None = visual.TextStim(win,text='None',color='black',pos=(letcol_xnone,choice_ynone),height=24)
    
    # Choosen color number per letter
    globals() ['choicecode_{}'.format(letters[0])] = visual.TextStim(win,color='black',pos=(letcoln_x1,choice_y1[0]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[1])] = visual.TextStim(win,color='black',pos=(letcoln_x1,choice_y1[1]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[2])] = visual.TextStim(win,color='black',pos=(letcoln_x1,choice_y1[2]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[3])] = visual.TextStim(win,color='black',pos=(letcoln_x1,choice_y1[3]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[4])] = visual.TextStim(win,color='black',pos=(letcoln_x1,choice_y1[4]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[5])] = visual.TextStim(win,color='black',pos=(letcoln_x1,choice_y1[5]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[6])] = visual.TextStim(win,color='black',pos=(letcoln_x1,choice_y1[6]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[7])] = visual.TextStim(win,color='black',pos=(letcoln_x2,choice_y2[0]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[8])] = visual.TextStim(win,color='black',pos=(letcoln_x2,choice_y2[1]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[9])] = visual.TextStim(win,color='black',pos=(letcoln_x2,choice_y2[2]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[10])] = visual.TextStim(win,color='black',pos=(letcoln_x2,choice_y2[3]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[11])] = visual.TextStim(win,color='black',pos=(letcoln_x2,choice_y2[4]),height=choice_col_size)
    globals() ['choicecode_{}'.format(letters[12])] = visual.TextStim(win,color='black',pos=(letcoln_x2,choice_y2[5]),height=choice_col_size)

    
    # INSTRUCTIONS
    instr.setText(welcome_txt)
    instr.draw()
    win.flip()
    event.waitKeys()
    
    # Show signs
    letter_sign.setAutoDraw(True)
    color_sign.setAutoDraw(True)
    choices_sign.setAutoDraw(True)
    
    # Start choosing until they are done
    no_final_choice = True
    trial=0
    while no_final_choice:
    #for t in range(trials):
        trial += 1   
        print('##### Trial {} #####'.format(trial))
        
        # Show list of stimuli:
        for l in letter_holder:
            l.setAutoDraw(True)
        for cc in colorcode_holder:
            cc.setAutoDraw(True)
        for cp in colorpatch_holder:
            cp.setAutoDraw(True)
    
        # Ask letter
        letter_ask.setAutoDraw(True)
        win.flip()            
        
        ## Record keypresses for letter display
        letter_resp = show_input(win,response_input=letter_input,type_input='letter',options=letters) # fix backspace
        letter_input.setAutoDraw(True)
        print('Letter response: ', letter_resp)
        
        # Ask color:
        color_ask.setAutoDraw(True)
        win.flip()
        
        ## Record keypresses for letter display
        color_resp = show_input(win,response_input=color_input,type_input='color',options=colorcode_input_str)
        print('Color response: ', color_resp)
        color_input.setAutoDraw(True)
        
        # Assign color to letter:
        pairs_dict[letter_resp] = color_resp
        
        # Show all choices choice at the right side (mix letter-color)
        for let in pairs_dict:
            if pairs_dict[let] is not '':
                col = int(pairs_dict[let])
                # SHow colored letter
                if (let in letters) and (col in colorcode_input):
                    #Colored letter
                    globals() ['choice_{}'.format(let)].setText(let)
                    globals() ['choice_{}'.format(let)].setColor(color_dict[col])
                    globals() ['choice_{}'.format(let)].setAutoDraw(True)
                    # Color number
                    globals() ['choicecode_{}'.format(let)].setText(col)
                    globals() ['choicecode_{}'.format(let)].setAutoDraw(True)
                else:
                    choice_None.setAutoDraw(True)

        win.flip()
        
        # Remove 'None' from dictionary
        no_let = []
        no_col = []
        for let, col in pairs_dict.items():
            if let not in letters:
                print('- {} is not a letter of interest.'.format(let))
                no_let.append(let)
            if (col not in colorcode_input_str) and (col is not ''):
                print('- {} is not a color of interest.'.format(col))
                pairs_dict[let] = ''
                print('- Now key {} that had non-color {} is empty again.'.format(let,col))
        
        for nl in range(len(no_let)):
            pairs_dict.pop(no_let[nl])
            print('- Non-letter {} removed from dictionary'.format(no_let[nl]))
        
        # Remove answers
        letter_ask.setAutoDraw(False)
        color_ask.setAutoDraw(False)
        color_input.setAutoDraw(False)
        letter_input.setAutoDraw(False)
        
        final_choice_ask.draw()
        win.flip()
        choice_None.setAutoDraw(False)
        
        keys=[]
        final_choice = event.waitKeys(keyList=['y','return'])
        print('- final choice:', final_choice[0])
        
        if final_choice[0] == 'y':
            
            #Select choosen colors and check empty ones:
            chosen_colors = []
            empty = 0
            for col in pairs_dict.values():
                chosen_colors.append(col)
                if col == '':
                    empty +=1
                    
            #Check repeated colors:
            seen = set()
            uniq = []
            duplic = []
            for c in chosen_colors:
                if c not in seen:
                    uniq.append(c)
                    seen.add(c)
                else:
                    duplic.append(c)
        
            # Finish only if there is no empty or repeated colors
            if len(duplic)>0 or empty>0:
                keep_choosing.draw()
                win.flip()
                core.wait(2)
            else:
                no_final_choice = False
    
    # Save data
    writer_count = 0  
    for letter in pairs_dict:      
        if pairs_dict[letter] is not '':
            colornum = int(pairs_dict[letter])
            colorRGB = color_dict[colornum]
            
            data_pairs.loc[writer_count] = [
                    subject_ID,             # subject
                    subject_yoked,          # subject yoked to
                    letter,                 # letter
                    colornum,               # Color code on the screen
                    int(colorRGB[0]),      #R in 255
                    int(colorRGB[1]),      #G in 255
                    int(colorRGB[2]),      #B in 255
                    ]  
            
            data_pairs.to_csv(output_filename,sep='\t')
            if subject_ID < 200:
                data_pairs.to_csv(output_filename_g2,sep='\t') # G2 gets G1's colors
            writer_count += 1
                
# Close-up
win.close()
core.quit()


