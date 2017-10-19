#!/usr/bin/python3

import pandas as pd
from Bio import SeqIO
import os

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
                print ("KeyError, key not found in All.faa, for %s\n"%(i))
                seqs.append("KeyError, key not found in All.faa")

        #handling no gene present 
        else:
            seqs.append("")

    # adding empty lines to merge list with dataframe successfully
    for i in range(4):
        seqs.append("")
    
    df["seq"] = seqs
    
    return df


def make_dir(outputdirname):
    """ making file if not already present
    """
    if not os.path.exists(outputdirname):
        os.makedirs(outputdirname)
    return

def print_out(df, fileName = "output.txt"):
    """ prints df to file with [fileName].csv
    """
    if not fileName.endswith(".csv"):
        fileName += ".csv" 
    
    df.to_csv(fileName, encoding = 'utf-8', sep = "\t")

    

    return fileName

def write_and_grab_files(fileName,df,grab=False):

    # write and grab file flag

    listFormattedDf= format_df(df)

    for pfam in listFormattedDf:
        filePath = fileName + pfam.name
        printedFile = print_out(pfam, filePath)
        print ("printed to %s.." % (filePath))        

        # grab file flag
        if grab == True:
            openDF = open_proteinortho_output(printedFile)
            dfSeq = get_seq(openDF)
            filePath = fileName + pfam.name
            print_out(dfSeq,filePath)
            print ("appended with seq %s.." % (filePath))
    return


def outputmanager(df,cutOff):
    #handles output functions.
    print ("%s families meeting criteria\n" % (df.shape[0]))

    ans = input("Initializing output manager..\n write pfam to file? -w:\n grab sequencesans and write pfam to file? -g\n exit? -e \n: " )


    # exit flag
    if ans == "e":
        raise SystemExit
    else:


        # case for making filenames from different kinds of settings
        # one cutoff setting
        if  type(cutOff[0]) != list:
            #making single folder name
            outFileName = "In_%s_Out_%s"%(cutOff[0],cutOff[1])
            make_dir(outFileName)
            path="./%s/"%(outFileName)

            if ans =="g":

                write_and_grab_files(path,df,True)

            elif ans == "w":

                write_and_grab_files(path,df,False)

        # multiple cutoff settings.
        elif  type(cutOff[0]) == list:

            #making a list of foldernames   
            outFileNames = ["In_%s_Out_%s"%(x[0],x[1]) for x in setting]
            # making file.
            for outFileName in outFileNames:
                make_dir(outFileName)
                path="./%s/"%(outFileName)

            # checks what mode is selected and
            # looping through all created filenames:
                if ans == "w" and ans =="g":

                    write_and_grab_files(path,df,True)

                elif ans == "w" and not ans =="g":

                    write_and_grab_files(path,df,False)

        else:
            print ("something went wrong with the settings..")
            raise SystemExit



