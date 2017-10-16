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

def format_df(df):
    #Transpose and split df based on pfam, returns multiple dfs.
  
    formattedFrames =[]
    ind = list(df.index)

    df = df.T

    for i in ind:
        
        dfC=df[i]
        formattedFrames.append(dfC)
    return formattedFrames

def print_out(df,fileName="output"):
    
    # prints df with fileName to .csv handles earlies .csv separately to get correct naming
    if fileName.endswith(".csv"):
        fileName=("Seq_"+fileName)
    else:
        fileName +=".csv" 
    
    try:
    	
    	df.to_csv(fileName, encoding='utf-8',sep="\t")
    	print ("printed to %s.."%(fileName))
    except:
    	print("could not print to %s"%(fileName))
    return fileName


def get_seq(df,allfasta="All.faa"):
    
    #gets sequence from All.fasta, based on the column of genes, handles misses and hits, not multiple hits.
    from Bio import SeqIO
    try:
    	record_dict = SeqIO.to_dict(SeqIO.parse(allfasta, "fasta"))
    except:
    	print("function get_seq in out.py: BioPython Could not read from %s"%allfasta)

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