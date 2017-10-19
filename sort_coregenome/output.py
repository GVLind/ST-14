#!/usr/bin/python3

import pandas as pd
from Bio import SeqIO

def open_proteinortho_output(myproject_fasta):    
    """	Argument: filename
		Purpose: convert [filename] to pandas dataframe 
		Returns: dataframe
	"""
  
    df = pd.read_csv(myproject_fasta, sep = "\t", names = ["strain","gene","seq"])
    
    return df

def format_df(df):
    """
    Transpose and split df based on pfam, returns multiple dfs.
    for output of all dfs to separate files
    """
    formattedFrames = []
    ind = list(df.index)
    df = df.T

    for i in ind:
        
        dfC = df[i]
        formattedFrames.append(dfC)

    return formattedFrames

def print_out(df, fileName = "output.txt"):
    """ prints df to file with [fileName].csv
    """
    if not fileName.endswith(".csv"):
        fileName += ".csv" 
    
    df.to_csv(fileName, encoding = 'utf-8', sep = "\t")

    print ("printed to %s.." % (fileName))

    return fileName

def get_seq(df,allfasta="All.faa"):
    """ gets sequence from All.fasta using SeqIO based on the column of genes
        in dataframe handles misses and hits, not multiple hits.
        offset by -4 fill with "" because no sequences for additiona
        data at the end.
    """
    record_dict = SeqIO.to_dict(SeqIO.parse(allfasta, "fasta"))
    
    dfGenes = df["gene"].tolist()
   
    seqs = []

    #offset because family information
    for i in dfGenes[:-4]:
        
        if i != "*":
            try:
                rec=record_dict[i]
                seqs.append(str(rec.seq))

            except KeyError:
                print ("KeyError, key not found in All.faa, printed as seq")
                seqs.append("KeyError, key not found in All.faa")

        #handling no gene present 
        else:
            seqs.append("")

    # adding empty lines to merge list with dataframe successfully
    for i in range(4):
        seqs.append("")
    
    df["seq"] = seqs
    
    return df


#def make_dir(outputdirname):
 #   if not os.path.exists(outputdirname):
  #      os.makedirs(outputdirname)

def outputmanager(df,cutOff):
    #handles output functions.
    ans = input("grab sequences? (y/e) : " )
    
    if ans != "e":

        listFormattedDf= format_df(df)
        
        
        for i in listFormattedDf:
            print (i)
            fileName = print_out(i,i.name)
            df = open_proteinortho_output(fileName)
            if ans == "y":
                dfSeq = get_seq(df)
                print_out(dfSeq,fileName)