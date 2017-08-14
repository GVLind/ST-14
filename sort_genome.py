#!/usr/bin/python

def FormatNames (string):

	tmpstring = string.replace(" ","")
	tmpstring = tmpstring.replace("(","")
	tmpstring = tmpstring.replace(")","")
	tmpstring = tmpstring.replace("/","")
	tmpstring = tmpstring.replace("[","")
	tmpstring = tmpstring.replace("]","")
	tmpstring = tmpstring.replace("+","")
	tmpstring = tmpstring.replace("\'","")
	tmpstring = tmpstring.replace("\"","")
	tmpstring = tmpstring.replace("\n","")
	outstring =tmpstring.replace(">","\n>")

	return outstring

def MakeProteinfamilyFiles (PFlist,header, Dbfile,inputpath):
	import subprocess

	subprocess.call(["mkdir","%s/%s"%(inputpath,header)])
	print ("grepping amino acids for proteinfamilies belonging to %s"%(header))


	for i in PFlist:
		print ("%s/%s"%(i[0],len(PFlist)), end="\r")

		meta = i[1][0:4]
		for j in i[2:]:
			if j != "*":

				p1 = subprocess.Popen(["grep","%s"%(j),"%s"%(Dbfile)],stdout=subprocess.PIPE)
				proteinfamilyfile = open("%s/%s/%s%s"%(inputpath,header,i[0],meta),"a")
				p2 = subprocess.Popen(["sed","s/###/\\n/g"],stdin = p1.stdout,stdout = proteinfamilyfile)
				
				proteinfamilyfile.close()
				p1.stdout.close()
	print ("%s files written to %s"%(len(PFlist),header))

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

tmpstring=""
with open (Dbfile,'r') as fi:

	for line in fi:
		if line.startswith(">"):
			line = line + "###"
		tmpstring=tmpstring + line
	tmpstring = FormatNames(tmpstring)

fi.close()

Dbfile = inputpath+"/NewDbfile.faa.tmp"

with open (Dbfile,"w+") as fo:
	fo.write(tmpstring)
fo.close

# creating list containing proteinfamilies

foldercontent = os.listdir(inputpath)
proteinfamilies = []
proteinfamilycounter = 0


for item in foldercontent:
	if item == "myproject.proteinortho":
		
		with open("%s/myproject.proteinortho"%(inputpath),"r") as proteinorthoOut:

			for line in proteinorthoOut:
				proteinfamilyWnumber=[]
				splitline=line.split('\t')
				newsplitline=[]



				for i in splitline:
					tmpstring = FormatNames(i)
					newsplitline.append(tmpstring)


				proteinfamily=newsplitline[3:]


				
				proteinfamilyWnumber.append("ProteinFamily_%s"%(proteinfamilycounter))
				
				proteinfamilycounter +=1

				for i in proteinfamily:
					i = i.replace("\n","")

					proteinfamilyWnumber.append(i)
				proteinfamilies.append(proteinfamilyWnumber)

		proteinorthoOut.close()

for i in proteinfamilies:
	stars = 0
	references = 0
	referenceName =[]
	strains = 0
	total_indicies=0
	for j in i:
		if j.startswith("*"):
			stars += 1
			total_indicies +=1
		elif j.startswith("ref"):
			references += 1
			total_indicies +=1
			format_j= j.split("_prot")[0]
			referenceName.append(format_j)

		elif j.startswith("S"):
			strains += 1
			total_indicies +=1

	metadata = [total_indicies,references,strains,stars,referenceName]

	i = i.insert(1,metadata)

#dividing proteinfamilies according to occurences in reference strains or not

#making a list of references.

Reflist=[]
for i in proteinfamilies[0][2:]:

	if isinstance(i,str) and not i.startswith("New"):
		i=i.replace("_formatted.faa","")
		i=i.replace(".faa","")
		i=i.replace("New","")
		Reflist.append(i)

#list all
#MakeProteinfamilyFiles(proteinfamilies,"All",Dbfile,inputpath)

#list all ST-14
ST14Pangenome = []
for i in proteinfamilies:
	if (i[1][2])>0:
		ST14Pangenome.append(i)
MakeProteinfamilyFiles(ST14Pangenome,"ST14Pangenome",Dbfile,inputpath)

#list containing reference
ContainingReference=[]
for i in proteinfamilies:
	if (i[1][1])>0:
		ContainingReference.append(i)
#MakeProteinfamilyFiles(ContainingReference,"Ref",Dbfile,inputpath)

#list notcontaining reference
NotContainingReference=[]
for i in proteinfamilies:
	if (i[1][1])==0:
		NotContainingReference.append(i)
#MakeProteinfamilyFiles(NotContainingReference,"No_Ref",Dbfile,inputpath)

#further dividing families:
#creating a list of where to find the reference strains.

AppendedReflist=[]
for i in Reflist:
	tmplist=[i]

	AppendedReflist.append(tmplist)
#print (AppendedReflist)

for List in proteinfamilies:
	for Strains in List[1][4]:
		#print ("------------------")
		for Reference in Reflist:
			if Reference in Strains:
				
				for Refindex in AppendedReflist:

					if Refindex[0] ==Reference:
						Refindex.append(List[0])

# pulling the data of reference strains from the proteinfamilies list.

for Referencelist in AppendedReflist:
	Temporary_proteinfamily =[]
	for Reference in Referencelist:
		#print (Reference[1:])
		for proteinfamily in proteinfamilies:
			
			if proteinfamily[0] == Reference:
				Temporary_proteinfamily.append(proteinfamily)

	#MakeProteinfamilyFiles(Temporary_proteinfamily,Referencelist[0],Dbfile,inputpath)
