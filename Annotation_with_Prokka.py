#!/Usr/bin/python

def get_fasta_from_spades(InputPath,OutputPath):
	import os
	from shutil import copy
	
	#initiating variale	
	Cwd=os.getcwd()
	ParentPath=''
	ChildPath=''
	OutputFilePaths=[]
	
	#trying to make sense of bad path formatting:
	os.chdir(OutputPath)
	outputpath = os.getcwd()
	os.chdir(Cwd)

	#Checking if path is ok
	if os.path.isdir(InputPath):
	
		# Changing cwd to path
		os.chdir(InputPath)
	
		# Saves path as variable
		path = os.getcwd()
	
		# listing content of first gen directories (strains)
		parentdircontent = os.listdir(os.getcwd())
		
		for entry in parentdircontent:
			
			#checking if content is folder..
			if os.path.isdir('%s/%s'%(path,entry)):
				
				# Save path and move to it
				ParentPath = ('%s/%s/'%(path,entry))
				os.chdir(ParentPath)
				
				#list content of second gen. directories (spades output)
				childdircontent = os.listdir(os.getcwd())
				
				
				for entry in childdircontent:
					
					# Checks if dir
					if os.path.isdir('%s%s/'%(ParentPath,entry)):
						
						#Saves path and moves to it
						ChildPath = ('%s%s/'%(ParentPath,entry))
						os.chdir(ChildPath)

						# listing content of thrid gen. directories
						childdircontent=os.listdir(os.getcwd())

						for entry in childdircontent:

							# looking for scaffolds.fasta
							if entry == ('scaffolds.fasta'):
	
								# Formatting filename
								TargetPath=('%s%s'%(ChildPath,entry))
								Targetname=TargetPath.split('/')
								Targetname = ('%s-%s'%(Targetname[-2],Targetname[-1]))
							
								# creating copy to outputdir
								copydestination = ('%s/%s'%(outputpath,Targetname))
								copy(TargetPath,copydestination)
							
								# Save path to copied folder appends to list.
								OutputFilePaths.append(copydestination)
			else:
				print('isnochilddir')	
	else:
		print('isnoparentdirdir')
	# returns cwd to startingposition
	os.chdir(Cwd)

	# Returns list of paths to copied files.
	return (OutputFilePaths)

def format_spades_fastafiles (ScaffoldPaths,OutputPath):

	import os
	
	# Initiating	
	ORFDict={}
	FeedTheList = False
	NodeNum=1

	
	#trying to make sense of bad path formatting:
	Cwd=os.getcwd()
	os.chdir(OutputPath)
	outputpath = os.getcwd()
	os.chdir(Cwd)

	for Scaffolds in ScaffoldPaths:
		
		# Formats filename from path.
		Filenamefmt = Scaffolds.split('-')
		Filenamefmt = Filenamefmt[1]
		Filenamefmt = Filenamefmt.split('_')
		FileID = ('Strain_%s_Run_%s'%(Filenamefmt[1],Filenamefmt[0][0]))

		#Opens inputfiles from paths and creates a new outputfile simulatneously.
		with open(Scaffolds,'r') as inputfile, open ('%s/%s'%(OutputPath,FileID),'w') as outputfile:
			
			# loops thorugh content of inputfile
			for line in inputfile:
		
				# If startswith >, looks for headers in input file and writes a formatted one at corresponding position in new file.
				if line.startswith('>'):
					
					# Formatting of new header from path.
					Strainfmt = Scaffolds.split('-')
					Strainfmt = Strainfmt[1]
					Strainfmt = Strainfmt.split('_')
					NodeID = ('>%s_Run_%s_Node_%i\n'%(Strainfmt[1],Strainfmt[0][0],NodeNum))
					
					# Unique counter updated
					NodeNum+=1		
				
					# writes to newfile		
					outputfile.write(NodeID)
	
				# else case handels all non header lines i.e. strings od aminoacid.
				else:
					outputfile.write(line)

	#returns output path same as args.Output but verified as working, can be used as argument for Run_prokka
	return(OutputPath)


def Run_Prokka(Path):
	import os
	
	#Variables
	content =[]

	# moves to folder
	os.chdir(Path)
	
	# Listing content in current folder (Path)
	content = os.listdir(os.getcwd())
	
	# looping through content
	for i in content:

		# Checking for "strain", hardcoded earlier in script.
		if i.startswith('Strain_S'):

			# Formats filenames to fit tags in prokka.
			tmpname =i
			tmpname = i.split('_')
			x=tmpname[1]+'_'+tmpname[3]

			# Formats string for execution
			execstring =  ('prokka --outdir %s --locustag %s --prefix %s --force --mincontiglen 200 --evalue 1e-08 ./%s'%(x,x,x,i))
		
			# Executes in terminal
			os.system(execstring)	

# Since outdir is located in cwd "find" is used to collect all .faa files.
def collect_and_clean():
	import os
	
	# Creates new folder
	os.system ('mkdir Prokka_Output_AA ')

	# moves all .faa files to the new folder
	os.system ('find -name *.faa -exec cp {} ./Prokka_Output_AA \;')

	# change dir to newfolder and list its content
	os.chdir ('Prokka_Output_AA')
	content = os.listdir(os.getcwd())

	# loops through content
	for i in content:

		# opens old and newfile
		with open(i,'r') as oldfile, open ('New%s'%(i),'w') as newfile:

                        # loops thorugh content of oldfile
                        for line in oldfile:

                                # If line startswith ">" - looks for headers in oldfile and writes a formatted one at corresponding position in newfile.
                                if line.startswith('>'):

                                        # Formatting of new header from path. removing blankspaces
                                        Headfmt = line.replace(' ','_')

                                        # writes to newfile             
                                        newfile.write(Headfmt)

                                # else case handels all non header lines i.e. strings of aminoacids.
                                else:
                                        newfile.write(line)
	# Remove all Old files.
	os.system ('rm S*')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="InputDirectory", required=True, help= "SPAdes output file location")
parser.add_argument("-o", dest="OutputDirectory", required=True, help= "Desired output directory")
parser.add_argument("-c", dest="Cutoff", help="Integer for percent")
args = parser.parse_args()
parser.parse_args()



# 0 take parent directory of SPAdes output as $Argument: Input_directory.

# 1 open directoryy  of strain - take directory name  as variable.

# 2 Open one of 2 subolder

# 3 Copy .fasta file to $Argument: Output_directory

# 4 Convert .fasta file:	1 Change name to $ Argument: Input_directory
#				2 Change header name to something shorter than 20 chars. ORF 1 ... ORF X ?

# 5 For each .fasta file in $Argument: Output_directory run prokka script:
	# $ prokka --outdir [.fasta file == Argument: Output_directory] --locustag .fasta file --prefix .fasta file

# 6 Collect and clean name

# 0,1,2,3,4.1 
ScaffoldPaths = get_fasta_from_spades (args.InputDirectory, args.OutputDirectory)

# 4.2
path=format_spades_fastafiles (ScaffoldPaths, args.OutputDirectory)

# 5
# PATH from format_spades_fasta
Run_Prokka (path)

# 6
collect_and_clean()

