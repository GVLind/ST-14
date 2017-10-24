#!/usr/bin/env python3
"""
This program takes a tab delimited output file from protein ortho
and through user options set critera for gene families that 
should be printed as .csv files. Graphical representation
of the raw data is possible.

dependncies:
-pandas
-matplotlib
-Bio
usage: refactored.py [proteinortho output]

options during the program:
- drop unwanted strains
- choose outgroup
- set cutoff - for multiple cutoffs, hardcoding is the only option.
- print gene families meeting criteria.
- grab sequences and print gene families meeting criteria.

"""

import sys
import sort_ui
import output
import graphics
import sort_data_manipulation as data

cutOff,inputCsv,listOutGroup,verbosity = sort_ui.inputmanager(sys.argv[1])

"""
Area for cutoffs: change cutOff = [cutOff] to supress user input and 
use hardcoded settings.

example for single cutoff:
cutOff = [[1,0]]

for multiple cutoffs
cutOff = [[1,0],[0.93,0],[0.86,0],[0.79,0],[0.73,0],[0.66,0],[0.59,0]]
"""

cutOff= [cutOff]

for i in cutOff:

	OutDf,statsCsv = data.data_manipulations(i,inputCsv,listOutGroup,verbosity)

	output.outputmanager(OutDf, i)

	graphics.graph(statsCsv,cutOff)
