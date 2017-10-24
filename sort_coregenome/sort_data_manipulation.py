#!/usr/bin/python3
"""
performes calculations on a sorted data set
by setting proter id, calculate number of gene families meeting criteria,
and dropping rows not meeting criteria
"""

import pandas as pd

def set_index_name(df,index="pfam"):
    """ Purpose: set specific index translates index to pfam index +1
        Returns: indexed dataframe
    """
    indexList = list(df.index)
    newIndex = []

    for i in indexList:
        newIndex.append("%s%s" % (index, i + 1))

    I = pd.Index(newIndex)
    indexed_df = df.set_index(I)
    
    return indexed_df

def drop_rows_not_meeting_criteria(df, cutOff):
    """ dropping rows not meeting criteria specified by
        set_cutoff in sort_ui module.
    """
    cutIn = cutOff[0]
    cutOut = cutOff[1]

    # handles the cutoff for ingroup files.
    if cutIn == 1:
        df = df[df["qIn"] == cutIn]
    else:
        df = df[df["qIn"] > cutIn]

    # handles the cutoff for the outgroup
    if cutOut == 0:
        df = df[df["qOut"] == cutOut]
    
    elif cutOut == 1:
        df = df[df["qOut"] == cutOut]

    else:
        df = df[df["qOut"] < cutOut]

    return df

def quota_row_full_csv(df, outGroup):
    """	Argument: Data frame, list of outgroup
		Action:
        1 creates an ingroup by subtracting the outgroup from raw data
        2 converts dataframe to 0 and 1
        3 calculates the sum of each row and calculates the quota of hits/total
        4 appends hits and quota at the end of original datafram
        5 returns the original dataframe with appended calculated values.
	"""
    # 1
    headers = list(df.columns.values)
    inGroup = [x for x in headers if x not in outGroup]

    # 2
    dfN = (df != "*").astype(int)
    
    # 3 calculates all values and sum of elements of each row of ingroup and outgroup
    inHits  = (dfN[inGroup].sum(axis = 1))
    inTot   = (dfN[inGroup].count(axis = 1))
    outHits = (dfN[outGroup].sum(axis = 1))
    outTot  = (dfN[outGroup].count(axis = 1))
    
    # 4
    if outGroup == []:
        #adds quota to end of dataframe
        df['hitsIn'] = inHits
        df['qIn'] = (inHits / inTot)
        df['hitsOut'] = 0
        df['qOut'] = 0

    else:
        df['hitsIn'] = inHits
        df['qIn'] = (inHits / inTot)
        df['hitsOut'] = outHits
        df['qOut'] = (outHits / outTot)

    # 5
    return df

def data_manipulations(cutOff,InputCsv,listOutGroup,verbosity=1):
    """ handles the manipulation on the input file based on user settings passed
        from sort_ui. returns a CSV only with proteinfamilies that meets the specified criteria
        and a CSV with all protein families regardless if they meet the criterias.

        1 creating proper indexing
        2 calculated the number of present genes for each row in dataframe.
        3 check if rows meet user specified criteria. Drop if not
        4 returns Meeting criteria CSV which is a dataframe of only the strains that
        meet users criteria. statsCsv which is a dataframe of all evaluated data, regardless
        if it meets criteria.
    """
    # 1
    IndexedCsv = set_index_name(InputCsv)    
    # 2
    statsCsv = quota_row_full_csv(IndexedCsv,listOutGroup)
    # 3
    MeetingCriteriaCsv = drop_rows_not_meeting_criteria(statsCsv,cutOff)

    # output verbosity    
    if verbosity == 1:
        print("cutoff: ",cutOff)
        print("outgroup: ",listOutGroup)
        print("\n-----------input\n")
        print(InputCsv)
        print("\n-----------quotas added")
        print(statsCsv)
        print("\n-----------dropped rows not meeting criteria")
        print(MeetingCriteriaCsv)

    # 4 
    return MeetingCriteriaCsv,statsCsv

