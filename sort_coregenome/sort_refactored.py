#!/usr/bin/env python3


import sys
import sort_ui
import output
import graphics
import sort_data_manipulation as data

cutOff,inputCsv,listOutGroup,verbosity = sort_ui.inputmanager(sys.argv[1])

cutOff = [[1,0],[0.93,0],[0.86,0],[0.79,0],[0.73,0],[0.66,0],[0.59,0]]

for i in cutOff:

	OutDf,statsCsv = data.data_manipulations(i,inputCsv,listOutGroup,verbosity)

	output.outputmanager(OutDf, i)






#graphics.graph(statsCsv)


## some isseues here, can't pass a loop of settings, needs to be fixed
## commented out some of the code in output.oy

#print (cutOff)
#cutOff = [[1,0], [0.95,0], [0.90,0], [0.85,0], [0.80,0]]
#for i in cutOff:



