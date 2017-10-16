#!/usr/bin/python3



def open_proteinortho_output(myproject_fasta):    
    """	Argument: filename
		Purpose: convert [filename] to either a plain dataframe,
		or opens it with fixed names if it ends with .csv 
		Returns: dataframe
	"""
    file = myproject_fasta
    import pandas as pd
    
    if file.endswith(".csv"):
        df = pd.read_csv(file,sep="\t",names = ["strain","gene","seq"])
    else:
        df = pd.read_csv(file,sep="\t")
    
    return df


def drop_unncessary_info(DataFrame,columns = ['# Species', 'Genes', "Alg.-Conn."]):
	"""	Argument: DataFrame
		Purpose: Drops columns '# Species', 'Genes', "Alg.-Conn."
		Output: DataFrame
	""" 

	DataFrame.drop(columns, inplace=True, axis=1)
	return DataFrame


def make_list_of_int_from_string(outString,coulmns):
    
    """	Argument: outString,coulmns
		Purpose: make a list from a "," delimited string and return a list of 
		corresponding elements in columns
		Returns: list of selected columns
	"""
    outString=outString.split(",")
    outInt =[]
    outList = []
    try:   
        for i in range(len(outString)):
            outInt.append(int(outString[i]))

        for i in outInt:
            outList.append(coulmns[i])
    except:
        None
    return outList

def choose_outgroup(rawInputCsv):
    
    # interprets user input and returns a list of a defined outgroup
    coulmns = list(rawInputCsv.columns.values)
    ok = "n"
    
    while 1:    
        
        for i in range(len(coulmns)):
            print (i, coulmns[i])
        outString = input("\nOutgroup? delimit , :")
        
        outList = make_list_of_int_from_string(outString,coulmns)
        
        
        for i in outList:
            print (i)
            
        ok = input("\nCorrect? (y/n/e): ")
        
        if ok == "y":
            return outList
        
        elif ok == "e":
            raise SystemExit




def choose_drops(rawInputCsv):
    
    # interprets user input and returns a list of strains to be dropped
    columns = list(rawInputCsv.columns.values)
    ok = "n"
    dropList=[]
    while 1:    
        print("\n")
        for i in range(len(columns)):
            print (i, columns[i])
        outString = input("\nDrops? delimit , : ")
        
        if len(outString) == 0:
            return dropList
        dropList = make_list_of_int_from_string(outString,columns)
        
       
        for i in dropList:
            print (i)
            
        ok = input("\nCorrect? (y/n/e): ")
        
        if ok == "y":
            print ("\ndropped",dropList,"\n" )
            return dropList
        
        elif ok == "e":
            raise SystemExit
    

def drop_pfam(df,dropList):
    #drops columns from df based on input list
    
     #print("inside drop_pfam")
     
     df.drop(dropList, axis=1, inplace=True)
     return df

def inputmanager(inputFile='myproject.proteinortho'):
    # gluing together functions that manages the iunput
    
    rawInputCsv = open_proteinortho_output(inputFile)
    inputCsv = drop_unncessary_info(rawInputCsv)
    #here
  
    dropList = choose_drops(inputCsv)
  	

  	# here the input Csv is pruned according to user specification with drops
    if len(dropList)==0:
        prunedCsv = inputCsv
    else:
        prunedCsv = drop_pfam(inputCsv,dropList)
    # outgroup defined
    listOutGroup = choose_outgroup(prunedCsv)
    
        
    cutOff = 0.95
    verbosity = 0


    return cutOff,prunedCsv,listOutGroup,verbosity