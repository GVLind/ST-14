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

def write_and_grab_files(Path,df,grab=False):
    """ writes file using print_out function
        for arguments: Path is expected to contain path ./Foldername/
        df is expected to be a pandas DataFrame
        grab is expected to be boolean.

        two cases:
        1 for each line in dataframe:print file containing
        only protein names and strains.
        2. relies on 1, only executed if garab == True due to lack of
        programming knowledge first the a df i opened from the file from 1
        just printed, then the sequence is found using get_seq
        finally new file with appended sequence overwrites the old one.
    """
    listFormattedDf= format_df(df)

    # just stranin and gene
    for pfam in listFormattedDf:
        filePath = Path + pfam.name
        printedFile = print_out(pfam, filePath)
        print ("printed to %s.." % (filePath))        

        # grab file flag strain, gene and seq
        if grab == True:
            openDF = open_proteinortho_output(printedFile)
            dfSeq = get_seq(openDF)
            filePath = Path + pfam.name
            print_out(dfSeq,filePath)
            print ("appended with seq %s.." % (filePath))
    return


def outputmanager(df,cutOff):
    """handles output functions.
        df is expected to be a pandas dataframe
        cutOff either a list containing cutoff for out group and in group
        or a list of the abouve cutoffs.
        ans can take one of four states:
        -e exit
        -w write only strain and gene to folder based on cutoff
        -g grabs sequence and performs the same task as -w
        - everything else: causes a crash

        1. first the function checks if cutoff is list or list of lists
        2. filenames are generated. either single or multiple.
        3. for each filename a file is created using the filename from 2
        5. path is created from filename from 2
        4. based on user input files are printed to to the path of the new
           folder
    """


    print ("%s families meeting criteria\n" % (df.shape[0]))

    ans = input("Initializing output manager..\n write pfam to file? -w:\n grab sequencesans and write pfam to file? -g\n exit? -e \n: " )


    # exit flag
    if ans == "e":
        raise SystemExit
    else:


        # case for making filenames from different kinds of settings
        # one cutoff setting

        #if  type(cutOff[0]) != list:

        
            #making single folder name
        outFileName = "In_%s_Out_%s" % (cutOff[0],cutOff[1])
        make_dir(outFileName)
        path = "./%s/" % (outFileName)

        if ans =="g":

            write_and_grab_files(path,df,True)

        elif ans == "w":

            write_and_grab_files(path,df,False)
"""

not sure is need this if i pass one setting at the time from
sort_refactored.
        # multiple cutoff settings.
        elif  type(cutOff[0]) == list:

            #making a list of foldernames   
            outFileNames = ["In_%s_Out_%s" % (x[0],x[1]) for x in setting]
            
            # making file.
            for outFileName in outFileNames:
                make_dir(outFileName)
                path = "./%s/" % (outFileName)

            # checks what mode is selected and
            # looping through all created filenames:
                if ans == "g":

                    write_and_grab_files(path,df,True)

                elif ans == "w":

                    write_and_grab_files(path,df,False)

        else:
            print ("something went wrong with the settings..")
            raise SystemExit
"""


