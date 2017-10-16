#!/usr/bin/env python3


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


def data_manipulations(cutOff,listOutGroup,InputCsv,verbosity=0):

    statsCsv = quota_row_full_csv(InputCsv,listOutGroup)

    print (statsCsv)
    
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

cutOff,listOutGroup,indexedCsv,verbosity = sort_ui.inputmanager(sys.argv[1])
OutDf,statsCsv=data_manipulations(cutOff,listOutGroup,indexedCsv,verbosity)
graphics.graph(statsCsv,cutOff)
output.outputmanager(OutDf)


