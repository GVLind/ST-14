#!/usr/bin/python3

def flatten_FASTA_strings(filename):	# open and read fasta-files return each all ORFs as a list, 1 ORF per line.

	with open (filename, 'r') as fasta_in:			# Opens fastafile as fasta_in
		fastain=fasta_in.read()						# Reads line by line		
		fasta_out=[]								# declaring outfile

		fastain=fastain.replace(">","@@@>")			# In order to split the input string and still keep > as deliminer for sequence header /
													# > is changed to @@@> temporarily, then split at @@@ generating a list consisting of header and seq /
		fastain = fastain.split("@@@")				# for all found ORFs
													
		for line in fastain:						# looping into list of lines
													# cheking if > is present
			if line.startswith(">"):				
				line =line.replace("\n","\t",1)		# replacing first occurece of \n with \t, to be able to remove all \n later on
			
			line=line.replace("\n","")				# removing all \n, leaving the element of the list as just one line containing header + sequence delimited by tab
			fasta_out.append(line)					# appending the formatted line to the output list. 

		return fasta_out 							

def format_myprojectProteinortho_to_proteinfamily_dict (filename):	# open and format proteinortho output in a dictionary, keys = proteinfamily, value = formatted proteinortho output

	with open (filename, 'r' ) as myproject_proteinorho:				# open specified input file

		body 		=	myproject_proteinorho.read().split("\n")[:-1]	# body declared, splits at \n and removes last empty entry since each line in myproject.proteinorthoends with \n
		header 		=  	body[0].split("\t")								# declaring header as the top row of the myproject.proteinortho
		fmtfile		=	[]												# declaring variables x3
		pf			=	{}
		pf_count	=	1


		for line in body[1:]:											# for each line in body except header, == row 0
			line=line.split("\t")										# formatting lines to a list separated by splitting at \t. line length == header length
			fmtfile = []												# declaring empty list as placeholder later in the creation of the dict.
			for i in range(len(header)):								# loop for each entry in header which is exactly equal to line length, matching element in line to element in header

				new_element= [header[i],line[i]]						# new element for dictionary consists of header element + line element
				fmtfile.append(new_element)								# adding new element to placeholder.

			pf["pfam:"+str(pf_count)]=fmtfile							# adding placeholder list to dictionary, formats counter to str for correct concatenation str + str
			pf_count+=1

	myproject_proteinorho.close()										# closing file.
	return pf															# returning dictionary wih .key() = pfam:x and .value() = hits of ortholog in each strain 	#open and read myproject.proteinortho return dict of protein families

def get_sequence(filename,filenameAllfasta):	# takes output from format_myprojectProteinortho_to_proteinfamily_dict and appends sequence from All.faa

	
	pfamdict = format_myprojectProteinortho_to_proteinfamily_dict(filename) # get dict fom function: format_myprojectProteinortho_to_proteinfamily_dict 
	
	# manually create an All.faa in the same dir as the myproject.proteinortho
	# $ cat *.faa > All.faa

	flatAll = flatten_FASTA_strings(filenameAllfasta)
	count =0
	pfamcount = 0

	for pfam_key in pfamdict.keys():										# for each proteinfamily key, pfam_key in proteinfamily dictionary, pfamdict:
		pfam_list = pfamdict[pfam_key]										# get list of values for each key.

		for pfam_list_element in pfam_list[3:]:								# get all elements from list of values, not including first three digits.
			if pfam_list_element[1] != "*":									# only include hits, =! *
				for ORF in flatAll:

					if ORF.count(pfam_list_element[1]) > 0:
						pfam_list_element.append(ORF)
						flatAll.remove(ORF)
						count+=1
						#print ("ORFS found: ", count,"\r",end="")
						#print (pfam_list_element[1], "\n", ORF)
			else:
				pfam_list_element.append("*")
		pfamcount +=1
		print("pfams found: ",pfamcount,"/",len(pfamdict.keys()),"\r",end="")
     
	return pfamdict						

def make_csv_from_proteinfamily_dict(pfam_dict):	#makes a tab-delimited .csv file from dictionary with the format name=key.csv, content = value[0] \t value[2]


	for key in pfam_dict:
		with open(key+".csv","w") as csv_file:
			for value in pfam_dict[key][3:]:			#cut away first three list element in list of each dictionary .value()
				csv_file.write((value[0] + "\t" + value[2])+"\n")
				



# working test for getsequence
fmt_dict =get_sequence("myproject.proteinortho","All.faa")

for i in fmt_dict:
	print(i)
	for j in fmt_dict[i][3:]:
		for k in j:
			print(k)


#test for make_csv_from_proteinfamily_dict

#make_csv_from_proteinfamily_dict(get_sequence("testmyproject.txt","All.faa"))