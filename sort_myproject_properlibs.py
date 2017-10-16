#!/usr/bin/env python3


def drop_unncessary_info(DataFrame,columns = ['# Species', 'Genes', "Alg.-Conn."]):
	"""	Argument: DataFrame
		Purpose: Drops columns '# Species', 'Genes', "Alg.-Conn."
		Output: DataFrame
	""" 

	DataFrame.drop(columns, inplace=True, axis=1)
	return DataFrame

    
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

def set_index_name(df,index="pfam"):
    """	Argument: DataFrame
		Purpose: set specific index
		Returns: indexed dataframe
	"""
    import pandas as pd
    indexList = list(df.index)
    newIndex=[]
    for i in indexList:
        newIndex.append("%s%s"%(index,i+1))
    I = pd.Index(newIndex)
    
    indexed_df=df.set_index(I)
    
    return indexed_df

def set_in_out_group(listOut):
    """	Argument:
		Output:
		Returns:
	"""
    listOutGroup = ["PittEE","KR494"]
    return listOutGroup

def drop_rows_not_meeting_criteria(df,cutOff):
    #dropping rows not meeting criteria specified by set_reference
    df = df[df["qIn"] > 0.95]
    df = df[df["qOut"] < (1-cutOff)]
    return df

def quota_row_full_csv(df,outGroup):
    """	Argument:
		Output:
		Returns:
	"""
    headers = list(df.columns.values)
    inGroup = [x for x in headers if x not in outGroup]

    dfN = (df!="*").astype(int)
    #calculates all values and sum of elements of each row of ingroup and outgroup
    inHits  = (dfN[inGroup].sum(axis=1))
    inTot   = (dfN[inGroup].count(axis=1))
    outHits = (dfN[outGroup].sum(axis=1))
    outTot  = (dfN[outGroup].count(axis=1))
    
    #adds quota to end of dataframe
    df['hitsIn'] = inHits
    df['qIn'] = (inHits/inTot)
    df['hitsOut'] = outHits
    df['qOut'] = (outHits/outTot)

    return df

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


def drop_pfam(df,dropList):
    #drops columns from df based on input list
    
    df.drop(dropList, axis=1, inplace=True)
    return df

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
    
<<<<<<< HEAD

=======
def drop_pfam(df,dropList):
    #drops columns from df based on input list
    
     #print("inside drop_pfam")
     print(dropList)
     df.drop(dropList, axis=1, inplace=True)
     return df
>>>>>>> e2b41fb8b889795a96d5159e2336feb4132fb078


def print_out(df,fileName="output"):
    
    # prints df with fileName to .csv handles earlies .csv separately to get correct naming
    if fileName.endswith(".csv"):
        fileName=("Seq_"+fileName)
    else:
        fileName +=".csv" 
    print ("printing to %s.."%(fileName))
    df.to_csv(fileName, encoding='utf-8',sep="\t")
    return fileName


def get_seq(df,allfasta="All.faa"):
    
    #gets sequence from All.fasta, based on the column of genes, handles misses and hits, not multiple hits.
    from Bio import SeqIO
    record_dict = SeqIO.to_dict(SeqIO.parse(allfasta, "fasta"))
    
    dfGenes = df["gene"].tolist()
    seqs=[]
    for i in dfGenes:
        
        if i != "*":
            try:
                rec=record_dict[i]
               # print(rec.name,"\t",rec.seq)
                seqs.append(str(rec.seq))
            except KeyError:
                print ("KeyError, key not found in All.faa, printed as seq")
                seqs.append("KeyError, key not found in All.faa")
                

        else:
            seqs.append("")
    
    df["seq"]=seqs
    
    return df


def format_df(df):
    #Transpose and split df based on pfam, returns multiple dfs.
  
    formattedFrames =[]
    ind = list(df.index)

    df = df.T

    for i in ind:
        
        dfC=df[i]
        formattedFrames.append(dfC)
    return formattedFrames

def graph(df,cutoff):
    import matplotlib.pyplot as plt
    import numpy as np
    #rng = np.random.RandomState(10)  # deterministic random data


    fig, axes = plt.subplots(nrows=3)
    fig.subplots_adjust(hspace=0.4)
    ax0, ax1, ax2= axes.flatten()

    
    ax0.hist(df["hitsIn"], bins=35,color="darkgreen")  # arguments are passed to np.histogram
    ax0.set_title("Histogram of gene distribution of Out-group")
   
    
    ax1.hist(df["hitsOut"], bins=10,color="teal")  # arguments are passed to np.histogram
    ax1.set_title("Histogram of gene distribution of out-group")
    

    hb = ax2.hexbin(x=df['hitsIn'], y=df['hitsOut'],bins='log',gridsize=30,cmap="YlGn")
    cb = fig.colorbar(hb)
    cb.set_label('log10(N)')
    plt.show()
    


def data_manipulations(cutOff,listOutGroup,InputCsv,verbosity=0):

    statsCsv = quota_row_full_csv(InputCsv,listOutGroup)

    print (statsCsv)
    graph(statsCsv,cutOff)
    # evaluates if meeting criteria drop other rowa
    MeetingCriteriaCsv = drop_rows_not_meeting_criteria(statsCsv,cutOff)

    # output verbosity    
    if verbosity ==1:
        print("cutoff: ",cutOff)
        print("outgroup: ",listOutGroup)
        print("\n-----------input\n")
        print(InputCsv)
        print("\n-----------numerical")
        print(numericalInputCsv)
        print("\n-----------quotas added")
        print(statsCsv)
        print("\n-----------dropped rows not meeting criteria")
        print(MeetingCriteriaCsv)
        print("\n-----------expressed as original genes")
        print(pfamMeetingCriteria)
        
    return MeetingCriteriaCsv

def inputmanager(inputFile='myproject.proteinortho'):
    # gluing together functions that manages the iunput
    

    rawInputCsv = open_proteinortho_output(inputFile)
    
    inputCsv = drop_unncessary_info(rawInputCsv)
    
    indexedCsv = set_index_name(inputCsv)
    
    dropList = choose_drops(indexedCsv)
    
    if len(dropList)==0:
        prunedCsv = indexedCsv
    else:
        prunedCsv = drop_pfam(indexedCsv,dropList)
    
    listOutGroup = choose_outgroup(prunedCsv) #     ["PittEE","KR494"]
        
    cutOff = 0.95
    
    df=data_manipulations(cutOff,listOutGroup,indexedCsv,verbosity = 0)
    
    return df

def outputmanager(df):
    #handles output functions.
    ans = input("grab sequences? (y?) : " )    
    listFormattedDf= format_df(df)
    
    
    for i in listFormattedDf:
        fileName = print_out(i,i.name)
        df = open_proteinortho_output(fileName)
        if ans == "y":
            dfSeq = get_seq(df)
            print_out(dfSeq,fileName)
    
    return

import sys
df = inputmanager(sys.argv[1])
outputmanager(df)


