#!/usr/bin/python3

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="Inputfolder", required=True, help="folder containging myproject.proteinortho and additional .faa used creating myproject.proteinortho")
args = parser.parse_args()
parser.parse_args()


# compressing .faa files to a temporary database.

inputpath = args.Inputfolder
Dbfile = inputpath+"/Db.faa.tmp"
execstring =  "cat %s/*.faa > %s"%(inputpath,Dbfile)

os.system(execstring)

#formatting Db

#Maybe some kind of buffer owerflow. hav to fix this
tmpstring=""
with open (Dbfile,'r') as fi:

	for line in fi:
		if line.startswith(">"):
			line = line + "###"
		tmpstring=tmpstring + line
	tmpstring = tmpstring.replace("-","")
	tmpstring = tmpstring.replace("(","")
	tmpstring = tmpstring.replace(")","")
	tmpstring = tmpstring.replace("/","")
	tmpstring = tmpstring.replace("[","")
	tmpstring = tmpstring.replace("]","")
	tmpstring = tmpstring.replace("+","")
	tmpstring = tmpstring.replace("\n","")
	tmpstring =tmpstring.replace(">","\n>")
	tmpstring =tmpstring.replace("###","\n>")
fi.close()

Dbfile = inputpath+"/NewDbfile.faa.tmp"

with open (Dbfile,"w+") as fo:
	fo.write(tmpstring)
fo.close

# creating dictionary containing proteinfamilies

foldercontent = os.listdir(inputpath)
proteinfamilies = {}
proteinfamilycounter = 0


for item in foldercontent:
	if item == "myproject.proteinortho":
		
		with open("%s/myproject.proteinortho"%(inputpath),"r") as proteinorthoOut:

			for line in proteinorthoOut:
				splitline=line.split('\t')
				proteinfamily=splitline[3:-1]
				proteinfamilycounter +=1
				proteinfamilies.update({proteinfamilycounter:proteinfamily})

proteinorthoOut.close()


# grabbing sequence data

# Not working corectly, f-ing cal f-ing up my .txtoutputs
with open(Dbfile,"r") as Database:
	with open ("tmp.txt","w")as tmp:

		#print proteinfamilies[12]

		for key in proteinfamilies.keys():
			
			for val in proteinfamilies[key]:
				
				print val[:20]
				for line in Database:
					#print val[:15]

					print line[1:20]
					#print value
					#tmp.write("%s %s"%(val,line))
					#print val,line
					if val[:3]==line[1:3]:
						print line
						print val

tmp.close
Database.close

#>S9_21_01729_Elongation_factor_Tu_2
#S17_21_01615_Phosphoribosylformylglycinamidine_synthase
#S17_21_01615_Ph