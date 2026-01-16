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
sp  = 15.0  # raise y of symbol position letters
lh2 = 120   # letter size RSA arial black

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

