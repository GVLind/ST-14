#!/usr/bin/python3
########################################################################################################################
# managing output

def flatten_FASTA_strings(filename):												# open and read fasta-files return each all ORFs as a list, 1 ORF per line.
    #use biopython IOparser    
					
# Slow
    def get_sequence(pfamdict,filenameAllfasta):										
        # takes output from format_myprojectProteinortho_to_proteinfamily_dict and appends sequence from All.faa
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
        						pfam_list_element.append(ORF.split("\t")[1])		# splitting orf at tab to just append the sequence.
        						flatAll.remove(ORF)
        			else:
        				pfam_list_element.append("*")								# if no ORF is found appends "*
        		pfamcount +=1
        		print ("pfams found: ",pfamcount,"/",len(pfamdict.keys()),"\r",end="")
        	return pfamdict						
# working
    def format_myprojectProteinortho_to_proteinfamily_dict (filename):					# open and format proteinortho output in a dictionary, keys = proteinfamily, value = formatted proteinortho output
        with open (filename, 'r' ) as myproject_proteinorho:				# open specified input file
            import pandas as pd
            from pandas import DataFrame, read_csv
            return pf														# returning dictionary wih .key() = pfam:x and .value() = hits of ortholog in each strain 
        																		
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
        # working
        def prompt_interface():
        
        	print("\noutputmanager options: \n\n1 print to .csv files\n2 get sequences for hits\n3 print you sequences meeting criteria to the prompt\n4 print current settings to prompt\n0 exit\n")
        	userInput = input("what do you want to do? ")
        	print("\n")
        	return userInput
        # working
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
        # working
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
        # working
        def check_if_all_True(trueList):
        	#checks if all strains meets criteria in pfam.
        	print ("truelist",trueList)
        	#print ("truelist element 1",trueList[1])
        	
        	for i in trueList:
        		if type(i) != list:
        			print ("truelist not list")
        			return False 
        			
        	if trueList == None:
        		return None
        	elif all(element[1] == True for element in trueList):
        		return True
        
        	if trueList[0] == True:
        		return True
        	else:
        		return False
        	
        	if trueList[0][1] == True:
        		return True
        	else:
        		return False
        	return False
        def make_dict_of_pfam_meeting_criteria(settings,inputAllPfams,ListOfEvaluatedPfams): # makes dict new dict of pfams meeting criteria, discards rest.
        	# takes data from sort_pfam_evaluate_if_proteinfamily_meets_criteria(setting,pfam) for one single strain
        	# settings= get_settings(sys.argv[1])
        	# pfams = format_myprojectProteinortho_to_proteinfamily_dict(sys.argv[2])
        	#checking input parameters
        	if type(settings) != list:
        		print(" arg settings in function make_dict_of_pfam_meeting_criteria needs to be of type list, returning none")
        		return
        	elif type(inputAllPfams) != dict:
        		print(" arg inputAllPfams in function make_dict_of_pfam_meeting_criteria needs to be of type dict, returning none")
        		return
        	elif type(ListOfEvaluatedPfams) != list:
        		print(" arg ListOfEvaluatedPfams in function make_dict_of_pfam_meeting_criteria needs to be of type list, returning none")
        		return
        	for i in ListOfEvaluatedPfams:
        		if type(i[-1]) != str:
        			print("arg ListOfEvaluatedPfams in make_dict_of_pfam_meeting_criteria must be type sting, returning none")
        			return
        	allPfams = inputAllPfams
        	pfamsMeetingCriteria ={}
        	#for i in ListOfEvaluatedPfams:	#debug
        		#print (i)	#debug
        	#filter out pfam that is within the limits of the setting parameters:
        	#since numer of settings determine how many lists there will be there fore loop over range
        	
        	for i in ListOfEvaluatedPfams:
        		print (i)
        
        	for i  in range(len(ListOfEvaluatedPfams)):
        		currentSetting 	= settings[i]
        		evaluatedPfams	= ListOfEvaluatedPfams[i]
        		for evaluatedPfam in evaluatedPfams[:-1]:
        			pfam = evaluatedPfam[0]
        			pfamsMeetingCriteria[pfam] = allPfams[pfam]
        			pfamsMeetingCriteria[pfam].append(evaluatedPfams[-1])
        	return pfamsMeetingCriteria
        # working
        def sort_pfam_comparedigits_equal_to(value1,value2):								# logical comparison == 
        	pfam_status = False
        	if value1 == value2:	#checking if values match out cutoff and percent out
        		pfam_status = True
        		#print(value1,"=",value2 , pfam_status)	#debug
        		return pfam_status
        	else:
        		#print(value1,"=",value2 , pfam_status)	#debug
        		return pfam_status						
        # working
        def sort_pfam_comparedigits_less_than_and_equal_to(value1,value2):					# logical comparison <=
        	pfam_status = False
        	if value1 <= value2:	#checking if values match out cutoff and percent out
        		pfam_status = True
        		#print(value1,"<",value2 , pfam_status)	#debug
        		return pfam_status
        	else:
        		#print(value1,"<",value2 , pfam_status)	#debug
        		return pfam_status	
        # working
        def sort_pfam_comparedigits_bigger_than_than_and_equal_to(value1,value2):			# logical comparison >=
        	pfam_status = False
        	if value1 >= value2:	#checking if values match out cutoff and percent out
        		pfam_status = True
        		#print(value1,">",value2 , pfam_status)	#debug
        		return pfam_status
        	else:
        		#print(value1,">",value2 , pfam_status)	#debug
        		return pfam_status
        # working
        def sort_pfam_compare_occurence_settings(evaluation_list):							# takes composite list as argument, performs interpretation of settings data, passed for logical comparisons downstream
        	#checking input parameters
        	if type(evaluation_list) != list or len(evaluation_list) != 8:
        		print ("evaluation_list faulty, needs to be type list and contain exactly 8 elements, returning none")
        		return
        	elif type(evaluation_list[1]) != list or len(evaluation_list[1]) != 7:
        		print ("settings element of composite list needs to be of type list and contain exactly 7 elements, returning none")
        		return
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
        # working
        def format_setting_to_string (setting):
        	#takes one setting must be list
        	#checking parameters
        	if type(setting) != list:
        		print("arg setting format_setting_to_string type must be list, returning none")
        		return
        	if len(setting) != 7:
        		 print("arg setting format_setting_to_string type must have len 7, returning none")
        		 return
        	listToString=""
        	formattedSetting=""
        	#checks if element 3 is list and converts is to string
        	if len(setting[3])>1 and type(setting[3]) == list:
        		for i in setting[3]:
        			listToString +=(i+ " ")
        		setting[3]=listToString
        	#converts list to string element by element.
        	for i in setting:
        		formattedSetting+=(i+ " ")
        	return formattedSetting
        # working 
        def sort_pfam_evaluate_if_proteinfamily_meets_criteria (settings,pfam):				# for each pfam key, data is generated from the value, this is passed together with the used settings as a composite list for comparison. 
        	# takes settings and pfam as argument,
        	# settings = get_settings(sys.argv[1])
        	# pfams = format_myprojectProteinortho_to_proteinfamily_dict(sys.argv[2])
        	status_output 	= []
        	for keys in pfam:													#for each key in dictionary the composite list is created and passed for comparison.
        		#print("key %s evaluated"%(keys) ) 			#debug
        		count_in 		= 0
        		count_in_hit 	= 0
        		count_out 		= 0
        		count_out_hit 	= 0
        		inList 	=[]
        		outList =[]
        		#starting from third element because parameters from proteinortho still left.
        		for elements in pfam[keys][3:]:			# counts the number of hits and misses in myproject.proteinortho creates a unique composite list.		
        			#created an inlist and outlist based on settings.
        			for i in settings[3]:
        				if i not in elements[0]:
        					if elements[0] not in inList:
        						inList.append(elements[0])
        
        				elif i in elements[0]:
        					
        					if elements[0] not in outList:
        						outList.append(elements[0])
        
        		#counts hits and total number of strains in out and inlist, makes calculations.
        		for i in inList:			
        			count_in +=1
        			if elements[1]!="*":
        				count_in_hit+=1
        		for i in outList:			
        			count_out +=1
        			if elements[1]=="*":
        				count_out_hit+=1									
        		percent_out=100*count_out_hit/count_out
        		percent_in=100*count_in_hit/count_in
        		CompositeList=[keys,settings,count_in_hit,count_in,percent_in,count_out_hit,count_out,percent_out] # Here the composite list is created 		 
        		status = sort_pfam_compare_occurence_settings(CompositeList)								# List is passed for comparison settings - calculated data			
        		#can be unquoted for debugging don't forget #
        		print ("""\nSettings:\t%s\npfam:\t%s\nOutstrain:\t%s\nNo in:\t%s\nNo in hits:\t%s\npercent in:\t%s\nNo out:\t%s\nNo out hits:\t%s\npercent out:\t%s\nstatus:\t%s
        		#"""%(settings,keys,settings[3],count_in,count_in_hit,percent_in,count_out,count_out_hit,percent_out,status))
        	
        	
        		status_output.append([keys,status])
        
        
        	return status_output #passes evaluation data as a list. [Key, value=True or False.]
        def make_sorted_proteinfamily_dict(settings,pfam):									# makes two cases multiple strains or simple strain as outgroup, makes calculations for each setting:
        
        	#must handle multiple reference strains, if multiple check if all pfams are true, else discard.
        
        	ListOfEvaluatedPfam =[]
        
        	#print (settings) #debug
        	#one loop for each setting.
        	for setting in settings:
        
        		ListOfPfamAllTrue =[]
        
        			# sort away pfam that dont meet criteria 		
        		IdAndValue=sort_pfam_evaluate_if_proteinfamily_meets_criteria(setting,pfam)
        
        			#filters away all false hits
        
        		if check_if_all_True(IdAndValue):
        			for i in IdAndValue:
        				print ("all",i,setting)
        		
        
        		#if check_if_all_True(IdAndValue):
        		#	for i in IdAndValue:
        		#		print ("True",i,setting)
        			# converts settings to string
        		settingString = format_setting_to_string(setting)
        			# appends settings to the end of the pfams
        		IdAndValue.append(settingString)	
        
        			#passing settings to list
        		ListOfEvaluatedPfam.append(IdAndValue)
        
        	#print (ListOfEvaluatedPfam)	#debug
        	#sort_pfam_evaluate_if_proteinfamily_meets_criteria(setting,pfam) pass data to make the output dictionary of pfams that meet criteria
        	pfamsMeetingCriteria = make_dict_of_pfam_meeting_criteria(settings,pfam,ListOfEvaluatedPfam)		
        	return pfamsMeetingCriteria
        def input_output():
        
        	import sys
        	settings= get_settings(sys.argv[1])
        	pfams = format_myprojectProteinortho_to_proteinfamily_dict(sys.argv[2])
        "	pfamsMeetingCriteria = make_sorted_proteinfamily_dict(settings,pfams)
        "	output_manager(pfams,pfamsMeetingCriteria,settings)
        
        
        """
        def unit_test_check_if_allTrue():
        	testlist=[]
        	testlist.append([["x1",True],["x2",False],["x3",True]])
        	testlist.append([["y1",True],["y2",True],["y3",True]])
        	testlist.append([["z1",False],["z2",False],["z3",False]])
        	testlist.append([["a",True]])
        	testlist.append(["b",False])
        	testlist.append([])
        
        	referencelist = [False,True,False,True,False,None]
        	outputlist =[]
        	for i in testlist:
        		print (i)
        
        	for i in testlist:
        		print (i)
        		print (check_if_all_True(i))
        		outputlist.append(check_if_all_True(i))
        	if referencelist==outputlist:
        		return True
        	else:
        		return False
        """


#def unit_tester():
#	print ("starting unittester..")

#	print ("checking function check_if_all_True")
#	if unit_test_check_if_allTrue():
#		print ("ok")
#	else: 
#		print("unit_test_check_if_allTrue faild miserably")


#unit_tester()
#input_output()
# stuff to do:
	# out output.csv in correct folder.
	# make output with corresponding settings
	# debug

# Split further into functions
# Settings: Json Jaml