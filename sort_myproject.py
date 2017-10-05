#!/usr/bin/python3
########################################################################################################################
# managing output

def flatten_FASTA_strings(filename):												# open and read fasta-files return each all ORFs as a list, 1 ORF per line.


	with open (filename, 'r') as fasta_in:				# Opens fastafile as fasta_in
		fastain=fasta_in.read().rstrip()				# Reads line by line		
		fasta_out=[]									# declaring outfile

		#print(fastain)					#debug
		fastain = fastain.split(">")[1:]				#split and remove first empty list element	

														# for all found ORFs
		for line in fastain:							# looping into list of lines
			line =line.replace("\n","\t",1).replace("\n","")				# replacing first occurece of \n with \t, to be able to remove all \n later on
			
														# removing all \n, leaving the element of the list as just one line containing header + sequence delimited by tab
			fasta_out.append(line)						# appending the formatted line to the output list.

		return fasta_out 							

def get_sequence(pfamdict,filenameAllfasta):										# takes output from format_myprojectProteinortho_to_proteinfamily_dict and appends sequence from All.faa

	
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

		print ("pfams found: ",pfamcount,"/",len(pfamdict.keys()),"\r",end="")
     
	return pfamdict						

def format_myprojectProteinortho_to_proteinfamily_dict (filename):					# open and format proteinortho output in a dictionary, keys = proteinfamily, value = formatted proteinortho output

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

def make_csv_from_proteinfamily_dict(pfam_dict):									# makes a tab-delimited .csv file from dictionary with the format name=key.csv, content = value[0] \t value[2]


	for key in pfam_dict:												
		filekey=key.replace(":","")										#special case for the filename

		with open("%s.csv"%(filekey),"w") as csv_file:
			for value in pfam_dict[key][3:]:							#cut away first three list element in list of each dictionary .value()
				strainData 	= value
				strainId	= strainData[0]
				geneId		= strainData[1]

				if len(strainData) == 3:								#if len ==3 sequence data is present along with strain id and protein id 
					strainSeq = strainData[2]
					csv_file.write("{0}\t{1}\t{2}\n".format(strainId,geneId,strainSeq))
					

				elif len(strainData) == 2:								# if len ==2 value is thought to contain strain id and proteinid
					csv_file.write("{0}\t{1}\n".format(strainId,geneId))
			print ("file %s.csv created"%(filekey))

def prompt_interface():

	print("\noutputmanager options: \n\n1 print to .csv files\n2 get sequences for hits\n3 print you sequences meeting criteria to the prompt\n4 print current settings to prompt\n0 exit\n")
	userInput = input("what do you want to do? ")
	print("\n")
	return userInput

def output_manager(pfamDictIn,pfamDictOut,settings):
	

	dictIn 	= pfamDictIn
	dictOut = pfamDictOut

	while 1:
		userInput=prompt_interface()

		if userInput=="1":
			make_csv_from_proteinfamily_dict(dictOut)

		elif userInput == "2":
			allFasta = input("provide the name of your all.fasta:")

			dictOut=get_sequence(dictOut,allFasta)

			print ("Seqeunces succesfully fetched for you Ids.\n Your dictionary is now updated","blue")

		elif userInput == "3":
			import pprint
			pp = pprint.PrettyPrinter(indent=4)
			pp.pprint(pfamDictOut)
		elif userInput == "4":
			import pprint
			pp = pprint.PrettyPrinter(indent=4)
			pp.pprint(settings)
		elif userInput == "0":
			return
		else:
			print(" Choose one ot the above setings: ")


#########################################################################################################################
# managing data

def get_settings (reffile):															# parse the settings or reference file.

	with open (reffile,"r") as referencefile:							# opens reffile
		
		outlist = []													# declaring local variables
		settinglist = []
		settingsOut=[]

		rawlist = referencefile.read().replace("\n","")					#removing all newlines.
		rawlist = rawlist.split(",")									# splits lat ,

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
		#print (outlist[0]) 	#debug


		for element in settinglist:
			element=element.lower()										# making all lowercase
			subelements=element.split("_")								# splittint to subelements

			if subelements[3] == "all":									#handling case of all strains in outgroup
				subelements[3] = outlist[0]

			elif subelements[3].isdigit():								# handling case of single strain of output.
				subelements[3] = outlist[int(subelements[3])]			# changing from number to actual outstrains, this also passes the 0-case of the outlist..
			else:
				print ("cannot find outgroup number, #  of strain or elements is needed")
			
			settingsOut.append(subelements)								# appending all subelements to new settingsfile.

		return settingsOut												#returning the settingslist on the form 1: <,>,= for out, 2: % for out, 3: out , 4: strains(s) in out 5<,>,= for in, 6 % for in  7 , in. 

def sort_pfam_comparedigits_equal_to(value1,value2):								# logical comparison == 


	pfam_status = False

	if value1 == value2:	#checking if values match out cutoff and percent out
		pfam_status = True

		#print(value1,"=",value2 , pfam_status)	#debug
		return pfam_status
	else:
		#print(value1,"=",value2 , pfam_status)	#debug
		return pfam_status						

def sort_pfam_comparedigits_less_than_and_equal_to(value1,value2):					# logical comparison <=


	pfam_status = False
	if value1 <= value2:	#checking if values match out cutoff and percent out
		pfam_status = True

		#print(value1,"<",value2 , pfam_status)	#debug
		return pfam_status
	else:
		#print(value1,"<",value2 , pfam_status)	#debug
		return pfam_status	

def sort_pfam_comparedigits_bigger_than_than_and_equal_to(value1,value2):			# logical comparison >=


	pfam_status = False
	if value1 >= value2:	#checking if values match out cutoff and percent out
		pfam_status = True

		#print(value1,">",value2 , pfam_status)	#debug
		return pfam_status
	else:
		#print(value1,">",value2 , pfam_status)	#debug
		return pfam_status

def sort_pfam_compare_occurence_settings(evaluation_list):							# takes composite list as argument, performs interpretation of settings data, passed for logical comparisons downstream
	
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


	#checking outgroup															#calculated perccent out 	settings: out cutoff 				





	if evaluation_list[1][0] == "=":

		pfam_status = sort_pfam_comparedigits_equal_to(							int(evaluation_list[7]), 	int(evaluation_list[1][1]))

	elif evaluation_list[1][0] == "<":

		pfam_status = sort_pfam_comparedigits_less_than_and_equal_to(			int(evaluation_list[7]), 	int(evaluation_list[1][1]))

	elif evaluation_list[1][0] == ">":		
		pfam_status = sort_pfam_comparedigits_bigger_than_than_and_equal_to(	int(evaluation_list[7]), 	int(evaluation_list[1][1]))



	#checking ingroup															#percent in 				settings: in cutoff 
	if evaluation_list[1][4]  == "=" and pfam_status == True:
		pfam_status = sort_pfam_comparedigits_equal_to(							int(evaluation_list[4]), 	int(evaluation_list[1][5]))

	elif evaluation_list[1][4] == "<" and pfam_status == True:
		pfam_status = sort_pfam_comparedigits_less_than_and_equal_to(			int(evaluation_list[4]), 	int(evaluation_list[1][5]))

	elif evaluation_list[1][4] == ">" and pfam_status == True:
		pfam_status = sort_pfam_comparedigits_bigger_than_than_and_equal_to(	int(evaluation_list[4]), 	int(evaluation_list[1][5]))

	#print(pfam_status)
	return pfam_status

def sort_pfam_evaluate_if_proteinfamily_meets_criteria (settings,pfam):				# for each pfam key, data is generated from the value, this is passed together with the used settings as a composite list for comparison. 
	
	
	status_output 	= []

	for keys in pfam:
													#for each key in dictionary the composite list is created and passed for comparison.
		#print("key %s evaluated"%(keys) ) 			#debug
		count_in 		= 0
		count_in_hit 	= 0
		count_out 		= 0
		count_out_hit 	= 0
		

		for elements in pfam[keys]:			# counts the number of hits and misses in myproject.proteinortho creates a unique composite list.
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
		CompositeList=[keys,settings,count_in_hit,count_in,percent_in,count_out_hit,count_out,percent_out] # Here the composite list is created 
		
  
		status = sort_pfam_compare_occurence_settings(CompositeList)								# List is passed for comparison settings - calculated data
		
		status_output.append([keys,status])

										
	"""
	can be unquoted for debugging don't forget #
		print (#
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
				%(settings,
					keys,
					settings[3],
					count_in,
					count_in_hit,
					percent_in,
					count_out,
					count_out_hit,
					percent_out,
					status	))
	"""
	return status_output #passes evaluation data as a list. [Key, value=True or False.]
	
def sort_pfam_evaluate_if_proteinfamily_meets_criteria_multiple(settings,pfam):		# placeholder for later use
	print ("inside multiple")
	

	
	return				

def make_dict_of_pfam_meeting_criteria(settings,inputAllPfams,ListOfEvaluatedPfams): # makes dict new dict of pfams meeting criteria, discards rest.

	allPfams = inputAllPfams

	pfamsMeetingCriteria ={}


	#filter out pfam that is within the limits of the setting parameters:
	#since numer of settings determine how many lists there will be there fore loop over range

	for i  in range(len(settings)):
		currentSetting 	= settings[i]
		evaluatedPfams	= ListOfEvaluatedPfams[i]

		for evaluatedPfam in evaluatedPfams:
			pfam = evaluatedPfam[0]
			meetsCriteria = evaluatedPfam[1]
					
			if meetsCriteria == True:

				pfamsMeetingCriteria[pfam] = allPfams[pfam]
				#print(currentSetting, evaluatedPfam, meetsCriteria)		#debug				
			else:
				pass
				#print(currentSetting, evaluatedPfam, meetsCriteria)		#debug

	return pfamsMeetingCriteria

def make_sorted_proteinfamily_dict(settings,pfam):									# makes two cases multiple strains or simple strain as outgroup, makes calculations for each setting:
	
	outputPfamDictMeetingCriteria ={}
	ListOfEvaluatedPfam =[]

	for setting in settings:
		#print ("\nsetting evaluated:\n%s\n"%(setting))				#debug
																										#settings list[] contains a list if it's multiple inputfiles, else a string is present.

		if type(setting[3])==str:																		# if comparison with 1 refstrain
			IdAndValue=sort_pfam_evaluate_if_proteinfamily_meets_criteria(setting,pfam)
			
			ListOfEvaluatedPfam.append(IdAndValue)

		elif type(setting[3]==list):																	#if comparison multiple refstrains
			sort_pfam_evaluate_if_proteinfamily_meets_criteria_multiple(setting,pfam)
		
	pfamsMeetingCriteria = make_dict_of_pfam_meeting_criteria(settings,pfam,ListOfEvaluatedPfam)		#pass data to make the output dictionary of

	return pfamsMeetingCriteria



#test make sorted proteinfamily_dict

# stuff to do:
	# multiple ref-comparisons
	# output.csv
	# out output.csv in correct folder.


import sys

settings= get_settings(sys.argv[1])

#print (settings)

pfams = format_myprojectProteinortho_to_proteinfamily_dict(sys.argv[2])
pfamsMeetingCriteria = make_sorted_proteinfamily_dict(settings,pfams)





#output_manager(pfams,pfamsMeetingCriteria,settings)



"""
#test flatten_FASTA_strings
line = (flatten_FASTA_strings("testfasta.fasta"))
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

"""
#test flatten_FASTA_strings
line = (flatten_FASTA_strings("testfasta.fasta"))

import pprint
pp = pprint.PrettyPrinter(indent=4)

pp.pprint(line)
"""