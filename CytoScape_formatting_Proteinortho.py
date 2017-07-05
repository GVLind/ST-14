#!/usr/bin/python3

# Formats output from proteinortho "myproject.proteinortho" to .csv tab delimited Proteifamily - gene - Strain
def Format_project_proteinortho_to_Cytoscape_csv(ProteinOrthoProject):

	import os

	Firstline = []
	FirstLineFormatted = False
	PfCounter = 0
	with  open('NewtworkForCytoscape.csv','w+') as newfile:
		for line in ProteinOrthoProject:
			PfCounter += 1
			
			splitline = line.split('\t')
			splitline.remove(splitline[1])
			splitline.remove(splitline[1])
			splitline[-1]=splitline[-1][:-1]	
			splitline.remove(splitline[0])
			# Formatting lines to just contain number of strains containing the protein and the gene name, firstline formatted accordingly.	
			if  FirstLineFormatted == False:
				Firstline=splitline
				FirstLineFormatted = True
				
				# Formatting Firstline
			# takes all other lines (proteinfamilies) as is.	
			else:
				for index in range(len(splitline)):
					newfile.write('Protein_Family_%i\t'%(PfCounter))

					if splitline[index]==('*'):
												
						newfile.write('\n')
					else: 
						newfile.write('%s\t'%( splitline[index]))
						newfile.write('%s\n'%(Firstline[index]))
	newfile.close()


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="proteinorthofile", type=argparse.FileType('r'), required=True, help="myproject.proteinortho ")
args = parser.parse_args()
parser.parse_args()

Format_project_proteinortho_to_Cytoscape_csv(args.proteinorthofile)
