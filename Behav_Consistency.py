#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Letter-color Consistency test
O.Colizoli 2020
Each letter of the alphabet in random order x 2
Color wheel opens at a randomized color on each trial (but does not turn)
Python 2..7
"""
# data saved in ~/LogFiles/sub-XXX

# Import necessary modules
import random
import numpy as np
import pandas as pd
import os, time  # for paths and data
from IPython import embed as shell
import Tkinter as tk # py27
from tkColorChooser import askcolor

# Get subject number via tkinter (command line doesn't work in PsychoPy)
subject_ID = []
session = []
## INPUT WINDOW
class GetInput():
    def __init__(self):
        self.root2 = tk.Tk()
        self.root2.title("Subject and Session")
        # always put in same location
        w = 400 # width for the Tk root
        h = 200 # height for the Tk root
        # get screen width and height
        ws = self.root2.winfo_screenwidth() # width of the screen
        hs = self.root2.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/6) - (w/6)
        y = (hs/6) - (h/6)
        self.root2.geometry('%dx%d+%d+%d' % (w, h, x, y))
        # Subject
        self.e = tk.Entry(self.root2)
        self.e.insert(0, 'Subject Number')
        self.e.pack()
        self.e.focus_set()
        # Session
        self.e2 = tk.Entry(self.root2)
        self.e2.insert(0, 'Session')
        self.e2.pack()
        self.e2.focus_set()
        
        txt='If each letter of the alphabet\
        \nwere to have a unique color,\
        \nwhat color would it have?\
        \n\nThere are no right or wrong answers.'
        # instructions
        self.instr = tk.Label(self.root2, bg='white', text=txt, font=("Helvetica", 14))
        self.instr.pack()
        
        b = tk.Button(self.root2,text='OK',command=self.get_input)
        b.pack(side='bottom')
        
        self.root2.mainloop()
        
    def get_input(self):
        subj_str = self.e.get() 
        sess_str = self.e2.get()
        subject_ID.append(subj_str)
        session.append(sess_str)
        self.root2.destroy()
        
## ASK INPUT
app = GetInput()   # subject and session
subject_ID = int(subject_ID[0])
session = int(session[0])

## Create LogFile folder cwd/LogFiles
cwd = os.getcwd()
logfile_dir = os.path.join(cwd,'LogFiles','sub-{}'.format(subject_ID)) 
if not os.path.isdir(logfile_dir):
    os.makedirs(logfile_dir)
timestr = time.strftime("%Y%m%d-%H%M%S") 
output_alphabet = os.path.join(logfile_dir,'sub-{}_sess-{}_task-consistency_events_{}.tsv'.format(subject_ID,session,timestr))

### CONSISTENCY TASK ###
alphabet = ['a','b','c','d','e','f','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
alphabet = ['a','b','c']

REPS = 2 # number of times to repeat whole alphabet

RGBS = [] # save output
L = '2'  # place holder 

class Test():
    def __init__(self):
        self.counter = 1
        self.root = tk.Tk()
        self.root.title("Subject {} Session {}".format(subject_ID, session))
        # always put in same location
        w = 400 # width for the Tk root
        h = 300 # height for the Tk root
        # get screen width and height
        ws = self.root.winfo_screenwidth() # width of the screen
        hs = self.root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/6) - (w/6)
        y = (hs/6) - (h/6)
        self.root.geometry('%dx%d+%d+%d' % (w, h, x, y))
        self.open1 = tk.Button(self.root, text='Pick a color:', command=self.pick_a_color, font=('Helvetica', '36'),padx=5, pady=5)
        self.open1.pack(fill=tk.X, expand=False)    
        self.letter = tk.Label(self.root, bg='white', text=L, font=("Helvetica", 90))
        self.letter.pack()
        self.root.mainloop()
        
    def quit(self):
        RGBS.append( [L ,self.RGB, self.HEX, abc] )
        self.root.destroy()
            
    def pick_a_color(self,):        
        # GET COLOR CHOOSER NOT OPEN ON TOP OF ROOT
        self.RGB,self.HEX = askcolor((random.randint(0,255), random.randint(0,255), random.randint(0,255)), parent=None, title='Pick a color: {}'.format(L) )
        self.letter.configure(fg = self.HEX)
        if self.counter:
            exit_button = tk.Button(self.root, text='FINISHED', command=self.quit, font=('Helvetica', '28'))
            exit_button.pack()
            self.counter = 0
        self.root.mainloop()

# MAIN LOOP        
abc = 1 # round
for R in np.arange(REPS):
    random.shuffle(alphabet) 
    # Open a new GUI per letter        
    for L in alphabet:      
        app = Test()
    abc+=1

####################################
## SAVE OUTPUT & determine conditions
print(RGBS)

# Save colors
DFS = pd.DataFrame(RGBS)
DFS.columns = ["letter","rgb","hex","choice"]
DFS['subject'] = np.repeat(subject_ID,len(DFS))
DFS['r']            = [c[0] for c in DFS['rgb']]
DFS['g']            = [c[1] for c in DFS['rgb']]
DFS['b']            = [c[2] for c in DFS['rgb']]
DFS.to_csv(output_alphabet) # save all alphabet/preferences for both groups (also in case it goes wrong)
print('consistency test - success!')


##### OUTPUT FIGURE WITH COLORS #####
# Sort and show letters x 2 side by side
del tk  # py27
del askcolor
import matplotlib.pyplot as plt # doesn't work together with tkinter
import seaborn as sns
fig = plt.figure(figsize=(10,5))

# Sort so the same letters go side by side for each choice
try:
    DFS.sort_values(by=['choice', 'letter'],inplace=True)
except:
    DFS = DFS.sort(['choice', 'letter'])

DFS.reset_index(inplace=True)
for i,A in enumerate(alphabet):
    ax = fig.add_subplot(6,5,i+1)
    ax.text(0.5, 0.5, DFS['letter'][i], color=DFS['hex'][i],fontsize=18)
    ax.text(0.25, 0.5, DFS['letter'][i+len(alphabet)], color=DFS['hex'][i+len(alphabet)],fontsize=18)
    ax.set_axis_off()    

sns.despine(offset=10, trim=True)
plt.tight_layout()
fig.savefig(os.path.join(cwd,'LogFiles','sub-{}'.format(subject_ID),'sub-{}_colors.pdf'.format(subject_ID)))
print('success: sub-{}_colors.pdf'.format(subject_ID))

    
    
    
    
    
