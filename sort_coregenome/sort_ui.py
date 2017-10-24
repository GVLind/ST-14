#!/usr/bin/python3
"""
starts with inputmanager that calls other functions in order to
1 open input file
2 drop unncessary info 
3 choose drops, strains not included
4 choose outgroup
5 set cutoff levels
6 hardcoded verbosity
7 return
- cutOff a list of two elements float
- prunedCsv a csv file without unecessary info and strains.
- listOutGroup a list of the out group
- verbosity 1 or 0

"""



def open_proteinortho_output(myproject_fasta):
    """	Argument: filename
		Purpose: convert [filename] to either a plain dataframe,
		or opens it with fixed names if it ends with .csv 
		Returns: dataframe
	"""
    file = myproject_fasta
    import pandas as pd
    
    if file.endswith(".csv"):
        df = pd.read_csv(file, sep = "\t",names = ["strain", "gene", "seq"])
    else:
        df = pd.read_csv(file, sep = "\t")
    
    return df

def drop_unncessary_info(DataFrame, columns = ['# Species', 'Genes', "Alg.-Conn."]):
    """	Argument: DataFrame
    Purpose: Drops columns '# Species', 'Genes', "Alg.-Conn."
    Output: DataFrame
    """
    DataFrame.drop(columns, inplace = True, axis = 1)

    return DataFrame

def make_list_of_int_from_string(outString, coulmns):
    
    """	Argument: outString,coulmns
		Purpose: make a list from a "," delimited string and return a list of 
		corresponding elements in columns
		Returns: list of selected columns
	"""
    outString = outString.split(",")
    outInt = []
    outList = []

    try:   
        for i in range(len(outString)):
            outInt.append(int(outString[i]))

        for i in outInt:
            outList.append(coulmns[i])
    except:
        print("Error: function make_list_of_int_from_string could not create a list ")
    return outList

def choose_outgroup(rawInputCsv):
    """
    Interprets user input and returns a list of a defined outgroup
    """
    coulmns = list(rawInputCsv.columns.values)
    userans = "n"
    
    # starts a loop with two exit conditions. exit and return list.
    while 1:

        # print available strain names to STDOUT.
        for i in range(len(coulmns)):
            print (i, coulmns[i])

        # takes user input
        outString = input("\nOutgroup? delimit , :")
        
        # transforms user input string into a list.
        outList = make_list_of_int_from_string(outString, coulmns)
        
        # double check with user.
        for i in outList:
            print (i)
            
        userans = input("\nCorrect? (y/n/e): ")
        
        # interprets user input
        if userans == "y":
            return outList
        elif userans == "e":
            raise SystemExit
        else:
            None

def choose_drops(rawInputCsv):
    """
    interprets user input and returns a list of strains to be dropped
    """

    columns = list(rawInputCsv.columns.values)
    ok = "n"
    dropList = []

    while 1:
        # print available strain names to STDOUT.    
        print("\n")
        for i in range(len(columns)):
            print (i, columns[i])

        # interprets user input, also handles 0 drops straight away.
        outString = input("\nDrops? delimit , : ")
        if len(outString) == 0:
            return dropList

        # transforms user input string into a list.
        dropList = make_list_of_int_from_string(outString, columns)
        
        # double check with user.
        for i in dropList:
            print (i)
            
        ok = input("\nCorrect? (y/n/e): ")
        
        # interprets user input
        if ok == "y":
            print ("\ndropped", dropList, "\n" )
            return dropList
        
        elif ok == "e":
            raise SystemExit
    
def drop_pfam(df, dropList):
    """drops columns from pandas dataframe based on input list
    """
     
     df.drop(dropList, axis=1, inplace=True)
     return df

def set_cutoff(cutIn = 0.95, cutOut = 0.05):
    """ lets the user state the cutoffs as a number 100-0 input in pos 1 and output as pos 2.
        Return standard settings ifexit.
    """
    statement = ""
    while statement != "e":

        # attempting to parse user input into a list if successful: pass list of rules.
        statement = input ("""state cutoff 0-100 for in and out-put\ndelim space, exit e\n:  """)
        try:
            rule = [int(x) for x in (statement.split())] 

        except:
            print ("\ninput must be integer\n")
            # makes the while loop start over since next if can't handle an empty list.
            rule = []

        # attempting to devide by 100 if successful: return settings.
        if len(rule) == 2 and rule[0] < 101 and rule[1] < 101 :

            try:

                cutIn = rule[0]/100
                cutOut = rule[1]/100

                return [cutIn,cutOut]

            except:

                print("\nError gict two number between 0-100\n")
        else:
         print ("\nwrong format\n")

    # if exit case
    print("returning standard settings in: 0.95 out: 0.05\nexiting..")
    return [cutIn,cutOut]

def inputmanager(inputFile = 'myproject.proteinortho'):
    """ gluing together functions that manages the iunput
        could probably be managed with argsparse.
    """
    rawInputCsv = open_proteinortho_output(inputFile)
    inputCsv = drop_unncessary_info(rawInputCsv)
    dropList = choose_drops(inputCsv)

  	# here the input Csv is pruned according to user specification with drops
    # handles no drops.
    if len(dropList)==0:
        prunedCsv = inputCsv
    else:
        prunedCsv = drop_pfam(inputCsv,dropList)

    # outgroup defined
    listOutGroup = choose_outgroup(prunedCsv)
            
    cutOff = set_cutoff()
    verbosity = 0

    return cutOff, prunedCsv, listOutGroup, verbosity
