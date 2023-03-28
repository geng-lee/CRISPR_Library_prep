#!/usr/bin/python3

''' 

'''


import pandas as pd
import os
import argparse

argparser = argparse.ArgumentParser(description = 'This Software is a tool to design Guiding RNA. \n This scrip allow for the filtering based on off-target activity prediction.')
argparser.add_argument('-I','--OffTarget', metavar = 'file', dest = 'OffTarget', type = str, required = True, help = 'Off Target file')
argparser.add_argument('-C','--CountTreshold', metavar = 'integer', dest = 'N', type = int, default=5, required = False, help = 'Off Target count threshold')
argparser.add_argument('-T','--Treshold', metavar = 'integer', dest = 'Thres', type = int, default=1.0, required = False, help = 'Threshold CFD to consider')
argparser.add_argument('-O','--out', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'Ouput file name')


if __name__ == '__main__':
    args = argparser.parse_args()
    OffTarget=pd.read_csv(args.OffTarget,sep='\t',usecols=['cfdOfftargetScore','seqId','guideId','guideSeq'], dtype={'cfdOfftargetScore':float,'seqId':str,'guideId':str,'guideSeq'},na_values='None')
    df=OffTarget.loc[OffTarget.cfdOfftargetScore>=args.Thres]
    df.index=df.seqId+df.guideId
    filtered=df.index.value_counts()>=args.N
    Excluded=filtered.guideSeq[filtered]
    pd.Series(Excluded).to_csv(args.Output,sep='\t',header=None,index=None)
