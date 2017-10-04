#!/usr/bin/python3

def flatten_FASTA_strings(filename):					# open and read fasta-files return each all ORFs as a list, 1 ORF per line.



	with open (filename, 'r') as fasta_in:				# Opens fastafile as fasta_in
		fastain=fasta_in.read().rstrip()				# Reads line by line		
		fasta_out=[]									# declaring outfile

		print(fastain)					#debug
		fastain = fastain.split(">")[1:]				#split and remove first empty list element	

														# for all found ORFs
		for line in fastain:							# looping into list of lines
			line =line.replace("\n","\t",1).replace("\n","")				# replacing first occurece of \n with \t, to be able to remove all \n later on
			
								# removing all \n, leaving the element of the list as just one line containing header + sequence delimited by tab
			fasta_out.append(line)						# appending the formatted line to the output list.

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

		return pf														# returning dictionary wih .key() = pfam:x and .value() = hits of ortholog in each strain 

																		#open and read myproject.proteinortho return dict of protein families

def get_sequence(filename,filenameAllfasta):						# takes output from format_myprojectProteinortho_to_proteinfamily_dict and appends sequence from All.faa

	
	pfamdict = format_myprojectProteinortho_to_proteinfamily_dict(filename) # get dict fom function: format_myprojectProteinortho_to_proteinfamily_dict 
	
	# manually create an All.faa in the same dir as the myproject.proteinortho
	# $ cat *.faa > All.faa

	flatAll = flatten_FASTA_strings(filenameAllfasta)
	pfamcount = 0

	for pfam_key in pfamdict.keys():										# for each proteinfamily key, pfam_key in proteinfamily dictionary, pfamdict:
		pfam_list = pfamdict[pfam_key]										# get list of values for each key.

		for pfam_list_element in pfam_list[3:]:								# get all elements from list of values, not including first three digits.
			if pfam_list_element[1] != "*":									# only include hits, =! *
				for ORF in flatAll:

					if ORF.count(pfam_list_element[1]) > 0:					# if ORF is found in list of proteins
						pfam_list_element.append(ORF)
						flatAll.remove(ORF)


			else:
				pfam_list_element.append("*")								# if no ORF is found appends "*

		pfamcount +=1
		print("pfams found: ",pfamcount,"/",len(pfamdict.keys()),"\r",end="")
     
	return pfamdict						

def make_csv_from_proteinfamily_dict(pfam_dict):					# makes a tab-delimited .csv file from dictionary with the format name=key.csv, content = value[0] \t value[2]


	for key in pfam_dict:
		with open(key+".csv","w") as csv_file:
			for value in pfam_dict[key][3:]:							#cut away first three list element in list of each dictionary .value()
				csv_file.write((value[0] + "\t" + value[2])+"\n")
				
def get_settings (reffile):												# parse the settings or reference file.
	with open (reffile,"r") as referencefile:							# opens reffile
		
		outlist = []													# declaring local variables
		settinglist = []
		settingsOut=[]

		rawlist = referencefile.read().replace("\n","")					# splits lat ,
		rawlist = rawlist.split(",")

		if rawlist[len(rawlist)-1]=='':									# case that removes empty list element if edgecase with extra ,
			rawlist=rawlist[:-1]



		outlist.append("0")												# appending index 0 as reverve for special case all, see below.

		for item in rawlist:											# creating a list of the outgroup where 0 is all and digits >0 = numbered refrences.
			if item.startswith("#"):
				item=item.split("_")									# splits away #numer_ before index
				outlist.append(item[1])									# appends to outlist

			elif item.startswith("@"):						
				settinglist.append(item[1:])							# apends all settings to list.

		outlist[0]=(outlist[1:len(outlist)])							# when all elements of rawlist is sorted. 0-case, is added to outlist == all outstrains.


		for element in settinglist:
			element=element.lower()										# making all lowercase
			subelements=element.split("_")								# splittint to subelements

			subelements[3] = outlist[int(subelements[3])]				# changing from number to actual outstrains, this also passes the 0-case of the outlist..
			
			settingsOut.append(subelements)								# appending all subelements to new settingsfile.

		return settingsOut												#returning the settingslist on the form 1: <,>,= for out, 2: % for out, 3: out , 4: strains(s) in out 5<,>,= for in, 6 % for in  7 , in. 

def make_sorted_proteinfamily_dict(settings,pfam):					# makes two cases multiple strains or simple strain as outgroup
	
	for setting in settings:
		
		if type(setting[3])==str:
			stateofeval=sort_pfam_evaluate_occurence(setting,pfam)
			
			print("state of evaluate: ",stateofeval	)
			
		elif type(setting[3]==list):
			sort_pfam_multiple(setting,pfam)

def sort_pfam_multiple(settings,pfam):

	

	
	return				

				
def sort_pfam_evaluate_occurence (settings,pfam):					# evaluates the quota of hits / occurences of out our ingroup 
	


	for keys in pfam:
		count_in 	= 0
		count_in_hit 	= 0
		count_out 	= 0
		count_out_hit 	= 0

		for elements in pfam[keys]:
			if elements[0].find(settings[3])==-1:
				count_in +=1
				if elements[1]!="*":
					count_in_hit+=1
			if elements[0].find(settings[3])>-1:
				count_out +=1
				if elements[1]=="*":
					count_out_hit+=1
		percent_out=100*count_out_hit/count_out
		percent_in=100*count_in_hit/count_in

		list_eval=[]
		list_eval=(keys,settings,count_in_hit,count_in,percent_in,count_out_hit,count_out,percent_out,)
		status = sort_pfam_compare_occurence_settings(list_eval)
									#passes evaluation data as a list. 
		status_output.append[keys,status]		

		print ("""
		Settings: 	%s	
		pfam:		%s
		Outstrain:	%s
		No in:		%s
		No in hits:	%s
		percent in:	%s
		No out:		%s
		No out hits:	%s
		percent out:	%s 
		status: 	%s
				"""%(settings,
					keys,
					settings[3],
					count_in,
					count_in_hit,
					percent_in,
					count_out,
					count_out_hit,
					percent_out,
					status	))
	return status_output

def sort_pfam_compare_occurence_settings(evaluation_list):
	
	#Legend	-evaluation_list:	#Legend Settings:

	#keys,			[0]			#Out operator	[0]
	#settings,		[1]			#Out cutoff		[1]
	#count_in_hit,	[2]			#Out 			[2]
	#count_in,		[3]			#Out strain 	[3]
	#percent_in,	[4]			#In operator 	[4]
	#count_out_hit,	[5]			#In cutoff 		[5]
	#count_out,		[6]			#In 			[6]
	#percent_out,	[7]

	pfam_status = False


	#checking outgroup								#perccent out 				settings: out cutoff 				
	if evaluation_list[1][0] == "=":

		pfam_status = sort_pfam_comparedigits_equal_to(				int(evaluation_list[7]), 	int(evaluation_list[1][1]))

	elif evaluation_list[1][0] == "<":

		pfam_status = sort_pfam_comparedigits_less_than_and_equal_to(		int(evaluation_list[7]), 	int(evaluation_list[1][1]))

	elif evaluation_list[1][0] == ">":		
		pfam_status = sort_pfam_comparedigits_bigger_than_than_and_equal_to(	int(evaluation_list[7]), 	int(evaluation_list[1][1]))



	#checking ingroup								#percent in 				settings: in cutoff 
	if evaluation_list[1][4]  == "=" and pfam_status == True:
		pfam_status = sort_pfam_comparedigits_equal_to(				int(evaluation_list[4]), 	int(evaluation_list[1][5]))

	elif evaluation_list[1][4] == "<" and pfam_status == True:
		pfam_status = sort_pfam_comparedigits_less_than_and_equal_to(		int(evaluation_list[4]), 	int(evaluation_list[1][5]))

	elif evaluation_list[1][4] == ">" and pfam_status == True:
		pfam_status = sort_pfam_comparedigits_bigger_than_than_and_equal_to(	int(evaluation_list[4]), 	int(evaluation_list[1][5]))

	print(pfam_status)
	return pfam_status

def sort_pfam_comparedigits_equal_to(value1,value2):


	pfam_status = False

	if value1 == value2:	#checking if values match out cutoff and percent out
		pfam_status = True

		print(value1,"=",value2 , pfam_status)	#debug
		return pfam_status
	else:
		print(value1,"=",value2 , pfam_status)	#debug
		return pfam_status
def sort_pfam_comparedigits_less_than_and_equal_to(value1,value2):


	pfam_status = False
	if value1 <= value2:	#checking if values match out cutoff and percent out
		pfam_status = True

		print(value1,"<",value2 , pfam_status)	#debug
		return pfam_status
	else:
		print(value1,"<",value2 , pfam_status)	#debug
		return pfam_status
def sort_pfam_comparedigits_bigger_than_than_and_equal_to(value1,value2):


	pfam_status = False
	if value1 >= value2:	#checking if values match out cutoff and percent out
		pfam_status = True

		print(value1,">",value2 , pfam_status)	#debug
		return pfam_status
	else:
		print(value1,">",value2 , pfam_status)	#debug
		return pfam_status


"""
#test make sorted proteinfamily_dict
settings= get_settings("ref.txt")
pfams = format_myprojectProteinortho_to_proteinfamily_dict("testmyproject.txt")
sorted_pfams = make_sorted_proteinfamily_dict(settings,pfams)

print (sorted_pfams)


"""
"""
#test for set_ref

reference = get_settings("ref.txt")

for i in reference:
	print (i)
"""

"""
# working test for getsequence
fmt_dict =get_sequence("testmyproject.txt","All.faa")

for i in fmt_dict:
	print(i)
	for j in fmt_dict[i][3:]:
		for k in j:
			print(k)

"""
#test for make_csv_from_proteinfamily_dict

"""
make_csv_from_proteinfamily_dict(get_sequence("testmyproject.txt","All.faa"))
"""
line = (flatten_FASTA_strings("testfasta.fasta"))

import pprint
pp = pprint.PrettyPrinter(indent=4)

pp.pprint(line)