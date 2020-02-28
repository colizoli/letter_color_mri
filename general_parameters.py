#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
PARAMETERS - Letters, counterbalancing, behavioral lab
General parameters and stimuli
"""
# Behavioral lab (DCCN cubicles measured 4 Feb 2020)
# visual angle calculator: https://elvers.us/perception/visualAngle/
scnWidth, scnHeight = (1920, 1080)
screen_width        = 53.5 # centimeters 
screen_dist         = 55.0
white               = [255,255,255] # background screen color
grey                = [128,128,128]
black               = [0,0,0] # background screen color

lh  = 100 # letter size
fh  = 50  # fixation cross height
ww = 1000 # wrap width of instructions text

# font styles
font_english    = "Tahoma"
font_chinese    = "Hanzi-Kaishu"
font_hebrew     = "HebrewUniversal"   #"PCSB Hebrew"
font_trained    = "Arial Black"

# COUNTERBALANCING STIMULI
# Group 1 subjects 100,101,102 etc.
# Group 2 subjects 200,201,202 etc.
# Counterbalancing will be based on odd/even: 
# Trained set: np.mod(x,2)
# Untrained set: not np.mod(x,2)

# Letter sets
letter_sets_lower = [
    ['b','c','d','e','g','n','o','p','q','r','w','x','y'], # 0 False 
    ['a','f','h','i','j','k','l','m','s','t','u','v','z']  # 1 True 
]
letter_sets_upper = [
    ['B','C','D','E','G','N','O','P','Q','R','W','X','Y'], # 0 False
    ['A','F','H','I','J','K','L','M','S','T','U','V','Z']  # 1 True
]

# lower case and 8 most frequent out of the 13 traiend
letter_sets_lower_8_IM = [
    ['c','d','e','g','n','o','r','w'], # 0 False SET
    ['a','h','i','l','m','s','t','u']  # 1 True SET
]

letter_sets_IM_PRACTICE = ['A','B','C','D','E','F','G','H']

trained_colors = [
    [0,163,228],
    [161,199,70],
    [183,61,160],
    [181,44,67],
    [16,114,86],
    [237,114,162],
    [58,175,75],
    [248,154,28],
    [109,57,142],
    [239,79,41],
    [49,60,163],
    [255,211,0],
    [9,181,172]
]

## ICONIC MEMORY TASK
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
cue_target_color    = black  #'black' # target cue in red, holder cue in black
cue_holder_color    = [200,200,200] #gp.grey #'grey'
background_color    = white #'white' # code??
instructions_color  = black #'grey'
letters_color       = black # all black

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
