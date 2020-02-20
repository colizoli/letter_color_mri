#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Letter-color Consistency test
O.Colizoli 2020
Each letter of the alphabet in random order x 2
Color wheel opens at a randomized color on each trial (but does not turn)
"""
# data saved in ~/LogFiles/sub-XXX

# Import necessary modules
import random
import numpy as np
import pandas as pd
import os, time  # for paths and data
from IPython import embed as shell
try:
    import Tkinter as tk
except:
    import tkinter as tk
from tkinter.colorchooser import askcolor

# Get subject number (command line)
try:
    subject_ID = raw_input("subject number: ")  # Python 2
except:
    subject_ID = input("subject number: ")
subject_ID = int(subject_ID)

## Create LogFile folder cwd/LogFiles
cwd = os.getcwd()
logfile_dir = os.path.join(cwd,'LogFiles','sub-{}'.format(subject_ID)) 
if not os.path.isdir(logfile_dir):
    os.makedirs(logfile_dir)
output_alphabet = os.path.join(logfile_dir,'sub-{}_consistency.csv'.format(subject_ID))

alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
REPS = 2 # number of times to repeat whole alphabet

RGBS = [] # save output
L = '2'  # place holder 

class Test():
    def __init__(self):
        self.root = tk.Tk()
        self.root.geometry('400x200')
        self.counter = 1
        self.open1 = tk.Button(self.root, text='Pick a color:', bg="blue", command=self.pick_a_color, font=('Helvetica', '36'),padx=5, pady=5)
        self.open1.pack(fill=tk.X, expand=False)    
        self.letter = tk.Label(self.root, text=L, font=("Helvetica", 90))
        self.letter.pack()
        self.root.mainloop()
        
    def quit(self):
        self.root.destroy()
            
    def pick_a_color(self,):        
        self.RGB,self.HEX = askcolor((random.randint(0,255), random.randint(0,255), random.randint(0,255)), parent=None, title='Pick a color: {}'.format(L) )
        RGBS.append( [L ,self.RGB, self.HEX] )
        self.letter.configure(fg = self.HEX)
        if self.counter:
            exit_button = tk.Button(self.root, bg="blue",text='FINISHED', command=self.quit, font=('Helvetica', '28'))
            exit_button.pack()
            self.counter = 0
        self.root.mainloop()
# MAIN LOOP       
for R in np.arange(REPS):
    random.shuffle(alphabet) 
    # Open a new GUI per letter        
    for L in alphabet:      
        app = Test()

####################################
## SAVE OUTPUT & determine conditions
print(RGBS)

# Save colors
DFS = pd.DataFrame(RGBS)
DFS.columns = ["letter","rgb","hex"]
DFS['subject'] = np.repeat(subject_ID,len(DFS))
DFS['r']            = [c[0] for c in DFS['rgb']]
DFS['g']            = [c[1] for c in DFS['rgb']]
DFS['b']            = [c[2] for c in DFS['rgb']]
DFS.to_csv(output_alphabet) # save all alphabet/preferences for both groups (also in case it goes wrong)
print('consistency test - success!')

    
    
    
    
    