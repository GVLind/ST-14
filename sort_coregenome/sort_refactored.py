#!/usr/bin/env python3


import sys
import sort_ui
import output
import graphics
import sort_data_manipulation as data

cutOff,inputCsv,listOutGroup,verbosity = sort_ui.inputmanager(sys.argv[1])


OutDf,statsCsv = data.data_manipulations(cutOff,inputCsv,listOutGroup,verbosity)
output.outputmanager(OutDf, i)

#graphics.graph(statsCsv)


## some isseues here, can't pass a loop of settings, needs to be fixed
## commented out some of the code in output.oy

#print (cutOff)
#cutOff = [[1,0], [0.95,0], [0.90,0], [0.85,0], [0.80,0]]
#for i in cutOff:



