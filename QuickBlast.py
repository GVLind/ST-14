#!/usr/bin/python3

# Here the program takes the inputfile from -i and chops every entry indicated by '>'up to individual files in order to feed them separately to BLAST.
def TmpforBLASTdb (QueryInfile,verbosity):

	import os

	f 				= QueryInfile
	SeqDict			= {}
	Seq 			= ''
	NewFileNames= []
	currentprotein 	=''
	newline			=''

	
	print '-------------------------------------------\nCreating temporary files'
	
	for line in f:

		if line.startswith('>'):
			for char in line:
				
				if char == '>':
					newline= newline + char
				if char == ' 'or char == '_':
					newline = newline +'_'
				elif char.isalnum():
					newline= newline + char
				else:
					newline = newline +''

			currentprotein=newline
			Seq=''
			newline=''
			

		elif line[0].isalpha():
			line = line.replace('\n','')
			line = line.replace(' ','')
			Seq = Seq+line
			SeqDict[currentprotein] = Seq

	for filename in SeqDict.keys():
		NewFileName = 'Query_'+filename[1:]+'.tmp.fasta'
		TmpFile = open(NewFileName,'w+')
		TmpFile.write(filename+'\n'+SeqDict[filename])
		NewFileNames.append(NewFileName)

	if verbosity==True:
		print NewFileNames
	return (NewFileNames)

# here the program takes the collected data from proteinortho as sample creates a local Blast database using the os.system method.
# Have some hardcoded sequences in the execString that must be changed from case to case.
def MkBlastDB (sample,verbosity):
	import os

	print '-------------------------------------------\nRunning makeblastdb'
	database = str(sample) + '_DataBase'
	Output = 'BlastDb_' +database
	execString = ("makeblastdb -in "+ str(sample) + " -dbtype prot -out " + Output)
	if verbosity == True:
		print '-----------\nNOTE: Printing settings: ' +execString 
	os.system(execString)

	return (Output)

# Program looks for files with the suffix '.tmp.fasta' in the working directory and runs them one after one against a given Blast database
# Have some hardcoded sequences in the execString that must be changed from case to case.
def RunQueryAgainstSampleDatabse(QueryInfile,sample,verbosity):
	import os

	Nodes = TmpforBLASTdb(QueryInfile,verbosity)
	database = MkBlastDB(sample,verbosity)

	print '-------------------------------------------\nRunning Blast\n-----------\nNOTE: Printing settings'
	for QueryFile in Nodes:
		execString = ('blastp -query ' + str(QueryFile)+ ' -db '+str(database)+' -evalue 1e-10 -num_descriptions 35 -num_alignments 5 -num_threads 4 -out ' + 'Blast_'+str(QueryFile) +'_VS_givensample.ncbi.blastp')
		if verbosity == True:
			print 'Settings : %s\n' %(execString)
		os.system(execString)

	if verbosity==False:
		print '-----------\nNOTE: Keeping tmpfiles '
		os.system('rm *.tmp.fasta')

#Sniffs out the matched sequences in the .blastp file.
#Handles without any hardcoing, I think.
def HitSniffer(verbosity):

	import os
	#lists every file in working directory think i should use psubprocess instead, but i seems troublesome.
	Path=os.popen('pwd').read()
	FolderList=os.listdir(Path[:-1])
	HitDict = {}
	
	if verbosity == True:
		print '-------------------------------------------\nRunning HitSniffer function\n' 
	
	# loops through all files in working directory
	for file in FolderList:		

		
		# opens only with correct ending
		if file[-12:] == '.ncbi.blastp': 	
			with open(file,'r') as f:	

				# Initates variables that have to be resetted with each opened file
				HitList = []				
				Queryname = ''
				FeedTheList = False

				# Starts reading file line by line using
					
				for line in f:

					#looks for name of query sequence, resetted with every new file opened. makes it possible to add items to HitList from this point with FeedTheList= True'
					#Limits to 50 Chars
					if line.startswith('Query='):
						if len(line)<=50:
							Queryname = line[:-1]
						else:
							Queryname = '%s'%(line[:50])

					#Starts reading
					if line.startswith('Length'):

						FeedTheList = True

					# breaks "for loop" when hitting the first alignment starting with ">" this have to do with how blast handles .fasta file so this one i pretty constant i think.
					if line.startswith('>') == False:

						#Handling the No hits found situation.
						if line.startswith ('*****'):
							HitList.append('***** No hits found *****')
							break

						#Checks if it is legal to read from given point.
						if FeedTheList == True and line.startswith(' '):

							# Cram together the text, could have been done earlier. but couldn't be bothered atm.
							line = line.replace(' ','')

							# A lengthy "if not" case were I sort out everything but the nodes of the actual sequences producing significant alignment.
							# Could have been done between "Sequences producing significant alignments:" and the start of the alignments but then i had do design some other soultion for the no hits case aswell.
							if line.startswith('Score')==False and line.startswith('Sequence')==False and line.startswith('\n')==False and line.startswith('Query')==False:

								# splits off redundant information.
								line = line.split('#')
								HitList.append(line[0])
					else:
						break

				# Gather all data so far into a dictionary just before resetting Hitlist and queryname.
				QuerynameNew = Queryname.split('=')

				#Sorts list
				Hitlist = HitList.sort()
				#Adds Hitlist to Dictionary with protein as key and hitlist as value
				HitDict[QuerynameNew[1]] = HitList


	f = open('SniffedBlastHits', 'w')
	for key in HitDict.keys():
		f.write('@%s' % (key))
		for value in HitDict[key]:
			f.write('%s\n'%(value))
		f.write('\n')

	if verbosity == True:
		print '-----------\nNOTE:Printing hits From Blast\n'
		for key in HitDict.keys():
			print 'query: %s' % (key)
			for value in HitDict[key]:
				print '%s'%(value)
			print '\n'

	return(HitDict)

# Not recommended for large datasets i.e. many query files.
#Takes crazy long time if many query Seq is present
def FindBlastHitsInProteinFamily (ProteinOrthoDir,verbosity):

	import os
	from subprocess import Popen,PIPE, STDOUT,call


	print '-------------------------------------------\nRunning FindBlastHitsInProteinFamily function\n'
	
	BlastHits = HitSniffer(verbosity)

	if verbosity == True:			
			print '-----------\nNOTE: prints HitsSniffer output from outside function \n%s' %(BlastHits)

	# create data fromSniffed blasthits if dewfault dictionary from Hitsniffer
	if type(BlastHits) == dict:
		if verbosity == True:
			print '-----------\nNOTE: HitSniffer handover dict successful \n prints out results from HitSniffer function'
		content=[]
		for i in BlastHits:
			content.append('@%s'%(i[:-1]))
			if verbosity == True:
				print i
			for j in BlastHits[i]:
				content.append(j)
				if verbosity == True:
					print j
		
	else:
		print 'Hitsniffer not correctly executed - returning'
		return


	# Change working directory to childfolder contatining Protein Families
	# Hardcoded.
	os.chdir (ProteinOrthoDir)

	#Verbosity
	if verbosity==True:
		print '-----------\nNOTE: Changed cwd to \"%s\"'%(os.getcwd())
		print '------------\nNOTE: runs execString in terminal, Hit = query, finding = terminal output\n'
	
	# Initiates variables for creating list of findings using grep tool for each of the blast hits among ProteinOrthofiles
	SearchHitList = []
	SearchHitHeader = ''

	#Cycles through list of folder files.
	for Hit in content:
		
		# Filters the name of the protein that we are looking for and saves it in the variable "Hit".
		if Hit.startswith("@")==True:
			SearchHitHeader = Hit

		# finds entris in the inputfile with BlastHits that we want to find in Proteinfamily files.	
		elif Hit.startswith("@") == False and Hit.startswith('*')==False:

			#declares string that is going to be executed
			execString = 'grep -l %s ProteinFamily*' % (Hit)

			#Executes execString using Popen
			proc = Popen(execString, shell= True, stdout=PIPE)
			finding = proc.communicate()[0]

			# Sorts out findings that correlate to actual Protein Family files.
			if finding.startswith('Prot')==True:
				
				#Saves Hitheader = the protein that we are looking for and finding = the output from grep-command that finds the ProteiFamily file were the protein is present.

				SearchHitList.append([finding[:-1],SearchHitHeader[2:]])

				# Verbosity
				if verbosity==True:
					print 'execstring: %s \nHit: %s \nfinding: %s \n' % (execString, Hit,finding)
	
	#Change back to parent directory
	os.chdir('..')

	d={}
	tmp=[]
	for i in SearchHitList:
		tmp.append(i[1])
	
	newtmp =list(set(tmp))
	#print newtmp

	for i in newtmp:
		d[i] =[]


	
	for i in d:
		newnewtmp =[]
		for j in SearchHitList:
			if i ==j[1]:
				newnewtmp.append(j[0])
				newnewtmp =list(set(newnewtmp))
		d[i] = newnewtmp

	if verbosity==True:
		for i in d:
			print '%s\nfound in: \n%s'%(i, d[i])

	f = open('FoundProteinFamilies', 'w')
	for i in d:
		f.write('%s\nfound in: \n%s\n'%(i, d[i]))

	# Returns HitDictionary
	if verbosity==True:
		print '-----------\nNOTE: Returning ProteinFamily Hits"%s'%(d)
	return(d)

# Takes Dictionary from HitSniffer creates a dictionary with strains as number only
def MakeHitDict(SniffedDict):
	print '-------------------------------------------\nRunning: MakeHitDict'
	ShortHitDict={}

	for pair in SniffedDict:

		tmp = SniffedDict[pair]
		tmp.sort()

		NodeNameShortList=[]
		NodeNameNumberList=[]
		
		for NodeName in tmp:
			NodeNameSplitted = NodeName.split('_')
			NodeNameShort = NodeNameSplitted[0]
			NodeNameShortList.append(NodeNameShort)

		for S_NodeName in NodeNameShortList:
			if S_NodeName != ('***** No hits found *****'):

				NumberNodeName = int(S_NodeName[1:])
				NodeNameNumberList.append(NumberNodeName)
			else:
				NodeNameNumberList=[]
			
		NodeNameNumberList.sort()
		ShortHitDict[pair] = NodeNameNumberList

	#print  '-----------\nNOTE: Returning Hit Dictionary"\n%s'%(ShortHitDict)
	return (ShortHitDict)
#Takes a Dictionary of key:list and creates a matrix that is outputted to terminal and .xls file.
def MakeTableFromDictionary(DictionaryIn):

	print '-------------------------------------------\nMakeTableFromDictionary'

	Dictionary = DictionaryIn


	#print Dictionary
	Rowlist =[]

	for i in Dictionary:
		Rowlist.extend(Dictionary[i])

	Rowlist = list(set(Rowlist))
	Proteins = []
	Matrix = [[0 for x in range(len(Rowlist))] for y in range(len(Dictionary))] 

	mIndex = 0
	for pair in Dictionary.items():
		Proteins.append(pair[0]);
		temp = list(set(pair[1]));

		for i in temp:
			Matrix[mIndex][Rowlist.index(i)] = 1
		mIndex += 1
	
	#verbosity
	#for i in range(len(Matrix)):
	#	print Proteins[i][:-2],Matrix[i]

	#Writes to .xls file
	F=open('Matrix.xls','w')
	#Writes first Row
	F.write('\t')
	for i in Rowlist:
		F.write ('S%s\t'%(i))

	#Writes Rows correlating to Proteins 
	F.write('\n')
	for i in range(len(Matrix)):
		F.write ('%s\t'%(Proteins[i]))
		for i in Matrix[i]:
			F.write('%s\t'%(i))

		F.write('\n')

#Moves output from QuickBlast to foldertree
def SortOutputToFolderTree():
	import os
	import datetime
	
	#Creates timestamp to keep track of the files if several runs are performed
	timestamp=str(datetime.datetime.now())
	timestamp=timestamp.replace(' ','_')

	#Creates Foldersystem with 'Run' as parent directory
	execstring = 'mkdir BlastFiles_%s;mkdir BlastDb_%s;mkdir Run_%s'%(timestamp,timestamp,timestamp)
	os.system(execstring)


	#Checks if temporary files are present, makes a directory and moves all .tmpfiles to the new directory
	for i in os.listdir(os.getcwd()):
		if i.endswith('.tmp.fasta'):
			os.system('mkdir Query_tmp_%s'%(timestamp))
			os.system('mv *.tmp.fasta ./Query_tmp_%s'%(timestamp))
			os.system('mv Query_tmp_%s ./Run_%s'%(timestamp,timestamp))
			break

	# Moves Blast files to the correct directory DB in one Blast output in another. 
	os.system('mv BlastDb_* ./BlastDb_%s'%(timestamp))
	os.system('mv Blast_Query_* ./BlastFiles_%s'%(timestamp))

	# Moves specific outputfiles to Run directory
	os.system('mv FoundProteinFamilies ./Run_%s'%(timestamp))
	os.system('mv Matrix.xls ./Run_%s'%(timestamp))
	os.system('mv SniffedBlastHits ./Run_%s'%(timestamp))

	#Moves the blast-directories to run directory
	os.system('mv Blast* ./Run_%s'%(timestamp))	

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="infile", type=argparse.FileType('r'), required=True, help="Input query sequences in .fasta format that is going to be run against the local BLAST databse")
parser.add_argument('-s',dest="sample",required=True,help="Samplefile, or the collected parsed data from ProteinOrtho using proteinOrtho2fasta.py by Chayan")
parser.add_argument('-d',dest='ProteinOrthoDir',metavar='PROTEINORTHODIR',required=True,help='Dir of Protein Families')
parser.add_argument('-v',dest='verbosity',action='store_true',default=False,help='Increases output to terminal from function.')
args = parser.parse_args()
parser.parse_args()



#Full:
RunQueryAgainstSampleDatabse(args.infile,args.sample,args.verbosity)
#FindBlastHitsInProteinFamily(args.ProteinOrthoDir,args.verbosity)
SniffedDictionary=HitSniffer(args.verbosity)
MakeHitDict(SniffedDictionary)
MakeTableFromDictionary(MakeHitDict(SniffedDictionary))
SortOutputToFolderTree()

print '-------------------------------------------\nFinished'
