#!/usr/bin/env python3


def set_index_name(df,index="pfam"):
    """ Argument: DataFrame
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



def data_manipulations(cutOff,InputCsv,listOutGroup,verbosity=0):

    IndexedCsv = set_index_name(InputCsv)    
    
    statsCsv = quota_row_full_csv(IndexedCsv,listOutGroup)
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
        
    return MeetingCriteriaCsv,statsCsv


import sys
import sort_ui
import output
import graphics

cutOff,inputCsv,listOutGroup,verbosity = sort_ui.inputmanager(sys.argv[1])
OutDf,statsCsv=data_manipulations(cutOff,inputCsv,listOutGroup,verbosity)
graphics.graph(statsCsv,cutOff)
output.outputmanager(OutDf)


