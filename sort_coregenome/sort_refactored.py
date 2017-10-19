#!/usr/bin/env python3


import sys
import sort_ui
import output
import graphics
import sort_data_manipulation as data

cutOff,inputCsv,listOutGroup,verbosity = sort_ui.inputmanager(sys.argv[1])
OutDf,statsCsv=data.data_manipulations(cutOff,inputCsv,listOutGroup,verbosity)
output.outputmanager(OutDf,cutOff)

#graphics.graph(statsCsv)



#cutOff = [[1,0],[0.95,0.5],[0.9,0.10],[0.85,0.15]]


