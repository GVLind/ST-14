#!/usr/bin/python3 


def open_input_reference_sequence (RefSeqFile): # Takes argument reference sequence
        RefSeqAsList = []

        # Read in each line to a list.
        for i in RefSeqFile:
                RefSeqAsList.append(i)

        #Return list
        return (RefSeqAsList)

def make_formatted_file_from_list_and_outfmt(inlist,outfmt): # takes list and makes a new list formatted with outfmt parameter.
	
	NewFileList = []
	for fileline in inlist:
		fileline=fileline.replace('\n','')
		if fileline.startswith('>'):
			fileline=fileline.replace(' ','_')
			fileline = '\n>ref-%s_prot-%s\n'%(outfmt,fileline[1:])
			NewFileList.append(fileline)
		else:

			NewFileList.append(fileline)
	return NewFileList

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="InputRefSeq", type=argparse.FileType('r'), required=True, help="Input RefSeq ")
parser.add_argument("-of", dest="Outputformat", required=True, help="Output format ")
args = parser.parse_args()
parser.parse_args()



#Takes a file as input, have to be .fasta format.
#The idea is to take a proteome file from a NCBI database and format the sequence headers to fit downstream applications
# 1 The filename is declared.
# 2 The Output forumat is declared.
# 3 The file name is altered to "Outputformat"
# 4 Each sequence header is altered to 'Ref:_lowercase(Outputformat)_originalHeader'
# 5 Written to newfile: 'outfmt.faa'
# Done..

#1
FileAsList = open_input_reference_sequence(args.InputRefSeq) 

#2
OutFmt = args.Outputformat

#3 #4 
NewFileAsList = make_formatted_file_from_list_and_outfmt(FileAsList,OutFmt)

# 5
with open('%s_formatted.faa'%(OutFmt),'w+') as f:
	for line in NewFileAsList:
		f.write(line)
f.close()
