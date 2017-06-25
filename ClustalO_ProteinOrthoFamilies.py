#!/usr/bin/python3



def Open_ProteinFamily (ProteinFamily): # Takes argument ProteinFamily
	NodeList = []
	
	# Read in each line to a list.
	for i in ProteinFamily:
		NodeList.append(i[:-1])

	#Return list
	return (NodeList)



def Collect_Node_Data (ListOfProteinFamily,ProteinDirectory):# Open folder of Proteomes takes argument ProteinFamilyList 
	import os
	
	Path = ProteinDirectory
	ORFDict={}
	Dirstring=''
	Header=''		
	FeedTheList=False

	#Move into ProteinDirectory.
	StartDir= os.getcwd()
	print (StartDir)
	os.chdir(Path)
	print ('after %s'%(os.getcwd()))
	
	#For each line in ProteinFamilylist:
	for node in ListOfProteinFamily:
		
		for dirfile in os.listdir(os.getcwd()):
			if dirfile.endswith('.faa'):
			
				with open(dirfile,'r') as f:
					for l in f:

						if l.startswith('>%s'%(node)):
							FeedTheList=True
							Dirstring=''
							Header=l[:-1]

							continue

						elif  FeedTheList == True:

							
							if l.startswith('>'):
								FeedTheList=False
								ORFDict[Header] = Dirstring,dirfile 
								
							elif FeedTheList == True:
								Dirstring+=l[:-1]
	os.system('cd %s'%(StartDir))	
	return ORFDict

############################################################################################
#Run each enty of multifasta in ClustalO




import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-pf", dest="ProteinFamily", type=argparse.FileType('r'), required=True, help="Protein File, conatining NODES.. ")
parser.add_argument("-pd", dest="ProteinDirectory", required=True, help="Path to all .faa files from Prodigal ")

args = parser.parse_args()
parser.parse_args()


ListOfNodes= Open_ProteinFamily(args.ProteinFamily)

for i in ListOfNodes:
	print (i)

OrfDir = Collect_Node_Data(ListOfNodes,args.ProteinDirectory)

for i in OrfDir:
	print ('\n')
	print (OrfDir[i][1])
	print (i)
	print (OrfDir[i][0])



