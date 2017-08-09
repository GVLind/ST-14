#!/usr/bin/python3


# Inputtextfile manipulations 1/3
def Store_File_Line_By_Line(InputFile): #Takes file as inputfile and converts is into a list. Have to be stored in order to be used several times.

	FileList = []
	for l in InputFile:
		FileList.append(l)
	return (FileList)

 # Inputtextfile manipulations 2/3
def Find_Matches (filelist): # Takes List of files as input, preferrably from Store_File_Line_By_Line.

	from operator import itemgetter

	lines = filelist
	Mismatch =[]
	Match =[]
	All_Lines=[]
	
	#Reads Files line by line
	for line in lines:
		line =line.split('\t')
		
		# Reads all lines except first one.
		if line[0].startswith('#')==False:

			#Check if protein isn't present in all strains.
			#if != number of strains --> Mismatch
			#else --> Match
			if len(line)-3 != len(line [0]):

				# Makes indicies 1,2,3 to int. for easier sorting.
				int_str_line=[]
				for s in line:
					if s.isdigit():
						int_str_line.append(int(s))
					else:
						int_str_line.append(s)

				Mismatch.append(int_str_line)
				All_Lines.append(int_str_line)
			else:
				# Makes indicies 1,2,3 to int. for easier sorting.
				int_str_line=[]
				for s in line:
					if s.isdigit():
						int_str_line.append(int(s))
					else:
						int_str_line.append(s)

				Match.append(int_str_line)
				All_Lines.append(int_str_line)


	#return int for indicecs1 to 3, str for the rest..
	Sorted_Mismatch = sorted(Mismatch,key=itemgetter(0))
	Sorted_Match = sorted(Match,key=itemgetter(0))
	Sorted_All_Lines = sorted(All_Lines,key=itemgetter(0))

	return(Sorted_Mismatch,Sorted_Match,Sorted_All_Lines)#Splits ProteinOrtho file in three lists, sorts and returns

# Inputtextfile manipulations 3/3
def MisMatchWithCutoff(ManipulatedInputList,numberpercent): #Takes list as input, preferrably the All list. 

	# Sorts out the variable for the cutoff value given by user.
	AboveCutOff=[]
	UnderCutOff=[]
	CutOff=numberpercent  
	Operator = 1.0

	#Checks if digit, else set cutoff to 0% (displays all)
	try:
		if CutOff.isdigit():
			CutOff=float(numberpercent)
	except:
		CutOff =0.0 

	if CutOff<=1:
		Operator = CutOff
	elif CutOff>1 and CutOff<=100:
		Operator = CutOff/100
	else:
		print 'Give a percent value 0-1 or 1-100, set to 100%'

	print 'set to %f'%(Operator)


	#imports sorted lists of the Manipulated list
	Lists = ManipulatedInputList

	for i in Lists:
		if Operator <= (float(i[0])/(len(i)-3)):
			AboveCutOff.append(i)
		else:
			UnderCutOff.append(i)

	#Returns Lists of above cutoff, udner cutoff and the cutoff value.
	return (AboveCutOff,UnderCutOff,CutOff)

# Output  manipulations 1/3
def Make_directory_In_Cwd(DirName,InCutOff):
	import os
	import datetime
	
	#Creates timestamp to keep track of the files if several runs are performed
	timestamp=str(datetime.datetime.now())
	timestamp=timestamp.replace(' ','_')
	timestamp=timestamp.replace(':','.')

	dirname =str(DirName)
	CutOff = str(InCutOff)
	

	#Creates Foldersystem with 'Run' as parent directory
	execstring = ('%s_%s')%(dirname,timestamp)
	os.system('mkdir %s'%(execstring))


	return (execstring) #Returns Name of creaated file.

# Output  manipulations 2/3
def Create_Files_from_List(DirName,List):
	import os

	print 'Create_Files_from_List ' + DirName

	Path = os.getcwd()+'/'+DirName+'/'
	Counter = 1

	for i in List:
		
		tmplist=[]
		for j in range(3,len(i)):
			
			if i[j].startswith('*') == False:
				tmplist.append(i[j])
		f = open('%sProteinFamily_%s_Hits:_%s'%(Path,Counter,len(tmplist)),'w+')

		for entry in tmplist:
			f.write(entry+'\n')
		Counter +=1
	os.system('cd ..')
	print '%s files written to %s'%(Counter,Path)
	return

# Output  manipulations 3/3
def Make_Dir_Of_List(InListOfMatch,InListOfMisMatch,CutOff):

	ListOfMatches = InListOfMatch
	ListOfMisMatches = InListOfMisMatch
	CutOffMatches=int(CutOff)
	CutOffMismatches=100-int(CutOff)	

	try:
		if type(ListOfMisMatches) == list and type(ListOfMatches) == list:
			print 'Succesful handover for %s' % (type(ListOfMisMatches))
	except:
		print 'input not list%s,%s'%(type(ListOfMisMatches),type(ListOfMatches))
		return
	
	MkdirStringMisMatches = Make_directory_In_Cwd('MisMatches-presence_under_%s'%(CutOff),CutOffMismatches)
	MkdirStringMatches = Make_directory_In_Cwd('Matches-presence_over_%s'%(CutOff),CutOffMatches)
	
	Create_Files_from_List(MkdirStringMisMatches,InListOfMisMatch)
	Create_Files_from_List(MkdirStringMatches,ListOfMatches)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-po", dest="ProteinOrtho", type=argparse.FileType('r'), required=True, help="ProteinOrtho Myproject.proteinortho")
parser.add_argument("-c", dest="Cutoff", help="Integer for percent")
args = parser.parse_args()
parser.parse_args()


ListOfFileentry = Store_File_Line_By_Line(args.ProteinOrtho)
CutOff=(args.Cutoff)

MatchList = Find_Matches(ListOfFileentry)[2]

ListOfMatches = MisMatchWithCutoff (MatchList,CutOff) [0]
ListOfMismatches = MisMatchWithCutoff (MatchList,CutOff) [1]

Make_Dir_Of_List(ListOfMatches,ListOfMismatches,CutOff)

