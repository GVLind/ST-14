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
	tmpstring = tmpstring.replace(" ","")
	tmpstring = tmpstring.replace("-","")
	tmpstring = tmpstring.replace("(","")
	tmpstring = tmpstring.replace(")","")
	tmpstring = tmpstring.replace("/","")
	tmpstring = tmpstring.replace("[","")
	tmpstring = tmpstring.replace("]","")
	tmpstring = tmpstring.replace("+","")
	tmpstring = tmpstring.replace("\n","")
	tmpstring =tmpstring.replace(">","\n>")

fi.close()

Dbfile = inputpath+"/NewDbfile.faa.tmp"

with open (Dbfile,"w+") as fo:
	fo.write(tmpstring)
fo.close

# creating dictionary containing proteinfamilies

foldercontent = os.listdir(inputpath)
proteinfamilies = []
proteinfamilycounter = 0


for item in foldercontent:
	if item == "myproject.proteinortho":
		
		with open("%s/myproject.proteinortho"%(inputpath),"r") as proteinorthoOut:

			for line in proteinorthoOut:
				proteinfamilyWnumber=[]
				splitline=line.split('\t')
				proteinfamily=splitline[3:]

				proteinfamilycounter +=1
				proteinfamilyWnumber.append("ProteinFamily_%s"%(proteinfamilycounter))
				
				for i in proteinfamily:
					i = i.replace("\n","")

					proteinfamilyWnumber.append(i)
				proteinfamilies.append(proteinfamilyWnumber)


				
				
proteinorthoOut.close()

Dblines=[]
with open (Dbfile,"r") as Db:

	for line in Db:
		Dblines.append(line)
Db.close()

for ProteinFamily in proteinfamilies:
	PFlist=[]
	for protein in ProteinFamily[1:]:
		
		for line in Dblines:

			if line.startswith(">"+protein):
				print "hit" + ProteinFamily[0]
				print protein
				print line
				PFlist.append(line)

		with open (inputpath+"/"+ProteinFamily[0],"w")as PF:
			for i in PFlist:
				i=i.replace("###","\n")
				PF.write(i)
		PF.close()


		





