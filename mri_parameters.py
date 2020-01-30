#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
PARAMETERS MRI (localizers & RSA tasks)
"""

# conda install pyqt
# pip install pyserial
import os
import pandas as pd
from psychopy import event
import general_parameters as gp

# BUTTON BOXES PUT ON USB MODE, CAN LISTEN AS NORMAL KEYBOARD
# http://www.psychopy.org/api/hardware/forp.html
# BB_SERIAL_PORT_BUTTONS = "COM1"
# from psychopy import hardware
# forp = hardware.forp.ButtonBox(serialPort=1, baudrate=19200)
# respond = forp,getEvents(returnRaw=False, asKeys=False, allowRepeats=False)

"""
General parameters and stimuli
"""

trigger_key = ['5']
#buttons = ['b','y','g','r','w']
buttons = ['1','2','3','4']
scanner = False # everything coming in as keyboard

if scanner: # when you are really scanning 
    ## Trigger from SerialPort: com3 SerialPort number, baudrate 115200, parity none, databits 8, stopbits 1, ASCIICODE 97
    BB_SERIAL_PORT_TRIGGER = "COM3"
    TRIGGER_CODE = 97
else: # wait for a '5' from the keyboard
    TRIGGER_CODE = trigger_key
    
def waitForTrigger(clock,): 
    '''  Need to listen on the serial port for the MRI start signal,
    Also/either wait for any key press ('q' will terminate). 
    https://groups.google.com/forum/#!topic/psychopy-users/TqxZYXx8aU4
    ''' 
    event.clearEvents(eventType='keyboard') # remove any keys waiting in the queue 
    while True: 
        if scanner == True: # set during the initial imports at the beginning of this module 
            SPORT = serial.Serial(BB_SERIAL_PORT_TRIGGER,115200) # OPEN PORT: port, baudrate, timeout (none)
            if (SPORT.inWaiting() > 0): # check for an MRI trigger signal on the serial port 
                print('MRI trigger received! {}'.format(SPORT.read())) 
                if (SPORT.read() == TRIGGER_CODE): # the MRI 'start' code, TRIGGER CODE
                    respond = [[TRIGGER_CODE,clock.getTime()]]
                    break # begin the experiment  
        else:
            respond = event.waitKeys(keyList=trigger_key+['q'], timeStamped=clock)
            # also check for a keyboard trigger
            if len(respond) > 0: 
                if respond == ['q']: core.quit() #  escape allows us to exit 
                break # else begin the experiment 
    return respond
           
# BOLD baseline timing (in seconds)
t_baseline  = 12   # 12 sec  

scnWidth, scnHeight = (1920,1080) # MRI BOLDscreen, # (1024,768) MRI stim computer, #(1280, 1024) # dummy scanner settings
screen_width        = 53.5 # centimeters
screen_dist         = 70.0

lh  = 100   # letter localizer size
lh2 = 120   # letter size RSA arial black
wh  = 75    # word size
fh  = 50    # fixation cross height
sp  = 15.0  # raise y of symbol position hebrew letters
ss  = 1.5   # symbol scaling
ww = 1000   # wrap width of instructions text

#### LOCALIZERS ####
# Number of trials/blocks
trials_loc      = 13    # number of trials per block
reps_loc        = 12    # number of 'blocks_VWFA'

# Localizer timing (RSA timing in own file)
t_word      = 0.75   
t_fix       = 0.25   
IBI         = 1     


"""
VWFA localizer (not using)
"""
# words           = pd.read_csv(os.path.join('stimuli','VWFA_stimuli.csv'))
# fonts_VWFA      = [gp.font_english,gp.font_hebrew]
# blocks_VWFA     = ['English','Hebrew']

"""
Letter localizer
"""
blocks_letter = ['Letter','Symbol']
fonts_letter = [gp.font_english,gp.font_hebrew]            
              
"""
Color localizer
"""
blocks_color = ['Color','Black']
symbols = ["A","B","C","D","E","F","G","H","I","J","K","L","M"] # present in Hebrew font

colors_RGB = [
    [163,0,228],
    [161,70,199],
    [61,183,160],
    [67,44,181],
    [86,114,16],
    [114,237,162],
    [175,58,75],
    [28,154,248],
    [57,109,142],
    [79,239,41],
    [60,49,163],
    [0,211,255],
    [172,181,9]
]



