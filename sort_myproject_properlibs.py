#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 14:39:14 2017

@author: gvl
"""



def flatten_fasta_seqIO():
    None
def drop_unncessary_info(DataFrame):
    #dropping unnecessary columns.
    columns = ['# Species', 'Genes', "Alg.-Conn."]
    DataFrame.drop(columns, inplace=True, axis=1)
    return DataFrame

    
def open_proteinortho_output(myproject_fasta):    
    #open tab delimited csv file.
    import pandas as pd
    file = myproject_fasta
    df = pd.read_csv(file,sep="\t")
    return df

def set_reference():
    #set reference
    cutOff=0.5
    
    return cutOff
    
def set_in_out_group(listOut):
    #set outgroup
    listOutGroup = ["PittEE","KR494"]
    return listOutGroup
    
def convert_to_numerical(df):
    dfN = (df!="*").astype(int)
    return dfN

def drop_rows_not_meeting_criteria(df,cutOff):
    #dropping rows not meeting criteria specified by set_reference
    df = df[df["qIn"] > cutOff]
    df = df[df["qOut"] < (1-cutOff)]
    return df

def quota_row_full_csv(df,outGroup):
    #defines ingroup based on outgropup
    headers = list(df.columns.values)
    inGroup = [x for x in headers if x not in outGroup]

    
    #calculates all values and sum of elements of each row of ingroup and outgroup
    inHits  = (df[inGroup].sum(axis=1))
    inTot   = (df[inGroup].count(axis=1))
    outHits = (df[outGroup].sum(axis=1))
    outTot  = (df[outGroup].count(axis=1))
    
    #adds quota to end of dataframe
    df['qIn'] = (inHits/inTot)
    df['qOut'] = (outHits/outTot)


    return df
def get_gene_names(dfPruned,dfRaw):
    rows = list(dfPruned.index)
    dfOut = dfRaw.iloc[rows]
        
    return dfOut
    
def setting_parser_with_json():
    None
    
def data_manipulations(cutOff,listOutGroup,InputCsv,verbosity=0):
    
    
    numericalInputCsv=convert_to_numerical(InputCsv)
    statsCsv = quota_row_full_csv(numericalInputCsv,listOutGroup)
    MeetingCriteriaCsv = drop_rows_not_meeting_criteria(statsCsv,cutOff)
    pfamMeetingCriteria = get_gene_names(MeetingCriteriaCsv,InputCsv)
    #get original sequence
    
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
        
    return pfamMeetingCriteria

def make_list_of_int_from_string(outString,coulmns):
    
    outString=outString.split(",")
    outInt =[]
    outList = []
       
    for i in range(len(outString)):
        outInt.append(int(outString[i]))

    for i in outInt:
        outList.append(coulmns[i])
    
    return outList

def choose_outgroup(rawInputCsv):
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

def get_ref(df,all_fasta="All.faa"):
    from Bio import SeqIO
    fasta_dict = SeqIO.index(all_fasta, "fasta")
    
    indexList = list(df.index)
    for i in indexList:
        print (fasta_dict[i])
    
    #for i in indexList:
     #   print fasta_dict[i]
    
    
def inputmanager():
    
    rawInputCsv = open_proteinortho_output('testmyproject.txt')
    inputCsv = drop_unncessary_info(rawInputCsv)
    listOutGroup = choose_outgroup(inputCsv)#["PittEE","KR494"]     #
    cutOff = set_reference()
    
    
    df=data_manipulations(cutOff,listOutGroup,rawInputCsv,verbosity = 0)
    
    return df

def print_out(df,fileName="output.csv"):
    print ("printing to %s.."%(fileName))
    df.to_csv(fileName, encoding='utf-8')

    
def outputmanager(df):
    
    listOfColumns=list(df.index)
    for i in listOfColumns:
        
        print(i,"\n",df.loc[i])
        get_ref(df.loc[i])
        #print_out(df.loc[i],"%s.csv"%(i))
    
    print ("\n-----------output")

    return
    
df = inputmanager()
outputmanager(df)



#print(df)
#write and outputmanager that gets sequences and prints to csv.
#outputmanager(df)