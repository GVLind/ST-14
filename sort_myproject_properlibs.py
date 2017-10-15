#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 14:39:14 2017

@author: gvl
"""

def drop_unncessary_info(DataFrame,columns = ['# Species', 'Genes', "Alg.-Conn."]):
    #dropping unnecessary columns.
    DataFrame.drop(columns, inplace=True, axis=1)
    return DataFrame

    
def open_proteinortho_output(myproject_fasta):    
    #open tab delimited csv file.
    file = myproject_fasta
    import pandas as pd
    
    if file.endswith(".csv"):
        df = pd.read_csv(file,sep="\t",names = ["strain","gene","seq"])
    else:
        df = pd.read_csv(file,sep="\t")
    
    return df

def set_index_name(df):
    # set index names as pfam + index
    import pandas as pd
    indexList = list(df.index)
    newIndex=[]
    for i in indexList:
        newIndex.append("pfam%s"%(i+1))
    I = pd.Index(newIndex)
    
    indexed_df=df.set_index(I)
    
    return indexed_df

def set_reference():
    #set reference
    cutOff=0.95
    
    return cutOff
    
def set_in_out_group(listOut):
    #set outgroup
    listOutGroup = ["PittEE","KR494"]
    return listOutGroup
    
def convert_to_numerical(df):
    #converts hits to numerical value 1 and misses 0
    dfN = (df!="*").astype(int)
    return dfN

def drop_rows_not_meeting_criteria(df,cutOff):
    #dropping rows not meeting criteria specified by set_reference
    df = df[df["qIn"] > 0.95]
    df = df[df["qOut"] < (1-cutOff)]
    return df

def quota_row_full_csv(df,outGroup):
    #defines ingroup based on outgropup
    headers = list(df.columns.values)
    inGroup = [x for x in headers if x not in outGroup]

    dfN = (df!="*").astype(int)
    #calculates all values and sum of elements of each row of ingroup and outgroup
    inHits  = (dfN[inGroup].sum(axis=1))
    inTot   = (dfN[inGroup].count(axis=1))
    outHits = (dfN[outGroup].sum(axis=1))
    outTot  = (dfN[outGroup].count(axis=1))
    
    #adds quota to end of dataframe
    df['qIn'] = (inHits/inTot)
    df['qOut'] = (outHits/outTot)


    return df
def get_gene_names(df,dfRaw):
    
    rows = list(df.index)
  
    dfOut = dfRaw.loc[rows]

    return dfOut
    



def make_list_of_int_from_string(outString,coulmns):
    
    #makes a list of an inputstring delimited by ","
    outString=outString.split(",")
    outInt =[]
    outList = []
       
    for i in range(len(outString)):
        outInt.append(int(outString[i]))

    for i in outInt:
        outList.append(coulmns[i])
    
    return outList

def choose_outgroup(rawInputCsv):
    
    # interprets user input and returns a list of a defined outgroup
    coulmns = list(rawInputCsv.columns.values)
    ok = "n"
    
    while 1:    
        
        for i in range(len(coulmns)):
            print (i, coulmns[i])
        outString = input("please declare you outgroup, delimited with , :")
        
        outList = make_list_of_int_from_string(outString,coulmns)
        
        
        for i in outList:
            print (i)
            
        ok = input("is this correct(y/n), exit for exiting..")
        
        if ok == "y":
            return outList
        
        elif ok == "exit":
            raise SystemExit

def choose_drops(rawInputCsv):
    
    # interprets user input and returns a list of strains to be dropped
    columns = list(rawInputCsv.columns.values)
    ok = "n"
    
    while 1:    
        
        for i in range(len(columns)):
            print (i, columns[i])
        outString = input("any drops? delimited with , :")
        
        dropList = make_list_of_int_from_string(outString,columns)
        
        
        for i in dropList:
            print (i)
            
        ok = input("is this correct(y/n), exit for exiting..")
        
        if ok == "y":
            print ("droplist",dropList )
            return dropList
        
        elif ok == "exit":
            raise SystemExit
    
def drop_pfam(df,dropList):
    #drops columns from df based on input list
    print("inside drop_pfam")
    print(dropList)
    df.drop(dropList, axis=1, inplace=True)
    return df



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
            rec=record_dict[i]
            print(rec.name,"\t",rec.seq)
            seqs.append(str(rec.seq))
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

def data_manipulations(cutOff,listOutGroup,InputCsv,verbosity=0):
    
    #convert to numerical
    numericalInputCsv=convert_to_numerical(InputCsv)
    print (numericalInputCsv)
    
    #get stats of each row ie pfam
    statsCsv = quota_row_full_csv(InputCsv,listOutGroup)
    print (statsCsv)
    # evaluates if meeting criteria drop other rowa
    MeetingCriteriaCsv = drop_rows_not_meeting_criteria(statsCsv,cutOff)
    
    # using original names
    #pfamMeetingCriteria = get_gene_names(MeetingCriteriaCsv,InputCsv)

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
    
    prunedCsv = drop_pfam(indexedCsv,choose_drops(indexedCsv))
    
    listOutGroup = choose_outgroup(prunedCsv) #     ["PittEE","KR494"]
        
    cutOff = set_reference()
    
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


