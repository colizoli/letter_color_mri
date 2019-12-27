#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
RSA Letter-Color Training, Colored Text
O.Colizoli 2019
Text Name: Arial Black, Font Size 12
Single letter
"""
import os, time # for paths and data
from docx import Document
from docx.shared import Pt
from docx.shared import RGBColor
from IPython import embed as shell
import numpy as np
import pandas as pd
import general_parameters as gp
import docx_search_loop as sl

# -------------
# Character formatting is applied at the Run level. 
# Examples include font typeface and size, bold, italic, and underline.
# A Run object has a read-only font property providing access to a Font object. 
# A run's Font object provides properties for getting and setting the character formatting for that run.
# -------------
# See: https://github.com/aekoch/docx_formatting_aekoch
# -------------
# First replace the entire documents font name and size, save it, then open again to run the search loop for the particular letters


inFileName = "test_book.docx"
midFileName = "middle_book.docx"
outFileName = "test_book_processed.docx"
subject_ID = int(raw_input("Subject ID: ")) 

#letters = ['a','b','c']
#RGBS = [(255,0,0),(0,255,0),(0,0,255)]

font_name = gp.font_trained
font_size = 12

## LETTER SET DEPENDENT ON SUBJECT NUMBER
trained_letters = gp.letter_sets_lower[np.mod(subject_ID,2)] # lower case, TRAINED set!!

# Get subject-specific letter-color mapping, letters always have same order, colors are permuted
DFS = pd.read_csv(os.path.join('stimuli','S{}.csv'.format(subject_ID)))
print(DFS)
subj_colors = [[DFS['r'][i],DFS['g'][i],DFS['b'][i]] for i in range(len(DFS))]
subj_letters = list(DFS['letter'])

letters = subj_letters
colors = subj_colors
# mappings = list(zip(letters,colors)) 


# First, replace the entire document's font name and size, save it
doc = Document(inFileName) # open an existing book
for paragraph in doc.paragraphs:
    for r in paragraph.runs:
        font = r.font
        font.name = font_name
        font.size = Pt(font_size)
        font.color.rgb = RGBColor(0, 0, 0) # set all to black 
doc.save(midFileName) # overwrite

# Second run the letter-by-letter search and replace loop
doc = Document(midFileName)
for paragraph in doc.paragraphs:
    for i,search in enumerate(letters):
        print(sl.find_occurances_in_paragraph(paragraph, search))
        format_func = lambda x:x.font.color.__setattr__('rgb', RGBColor(colors[i][0],colors[i][1],colors[i][2]))
        for start in sl.find_occurances_in_paragraph(paragraph, search):
            sl.apply_format_to_range(paragraph, start, start + len(search), format_func) 
doc.save(outFileName)