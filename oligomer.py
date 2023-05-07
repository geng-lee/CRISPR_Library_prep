#!/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os
import argparse
import re
import pysam

argparser = argparse.ArgumentParser(
    description='This Software is part of a tool to design Guiding RNA. \n This script helps design oligomer/concatemere for array Crispr experiment ')
argparser.add_argument('-S', '--Study', metavar='name',
                       dest='Guide', type=str, required=True, help='Guide Sequence')
argparser.add_argument('-P', '--Positive', metavar='file',
                       dest='Vep', type=str, required=True, nargs='*', help='Variant effect predictor associated with Positive control library')
argparser.add_argument('-N', '--Negative', metavar='name',
                       dest='Negative', type=str, required=True, nargs='*', help='Negative controls ')
argparser.add_argument('-O', '--out', metavar='file', dest='Output',
                       type=str, required=False, default='out', help='Ouput file name')
argparser.add_argument('--BSMBI', metavar='file', dest='BSMBI',
                       type=str, required=False, default='CGTCTC', help='Restriction Site')
argparser.add_argument('-p', '--PRIMER_FILE', metavar='file',
                       dest='PRIMER', type=str, required=False, help='Primer File')
argparser.add_argument('-f', '--forward_primer', metavar=int, dest='Primer_pair', type=int,
                       required=False, default=50, help='Designation of the primer pair')
argparser.add_argument('-F', '--fragments', metavar='Annotation', dest='fragments', type=str, required=False, nargs='*', default=[
                       'CGTCTCACACCG', 'GTTTTGAGACGgactgcCGTCTCcCACCG', 'GTTTaGAGACGggactaCGTCTCgCACCG', 'GTTTcGAGACGcttctcCGTCTCtCACCG', 'GTTTgGAGACG'], help='')


def Assert_Fragments_BSMBI(BSMBI, fragments):
    flags = ''
    BSMBI_reverse = str(Seq(BSMBI).reverse_complement())
    if fragments[0][0:len(BSMBI)] != BSMBI:
        flags = flags + \
            'Fragment 1 doesn\'t contain restriction site (BSMBI) at its beginning \n'
    for i in range(1, len(fragments)-1):
        if fragments[i].find(BSMBI) == -1:
            flags = flags + 'Fragment  ' + \
                (i+1) + ' doesn\'t contain forward restriction site (BSMBI) \n'
        if fragments[i].find(BSMBI_reverse) == -1:
            flags = flags + 'Fragment  ' + \
                (i+1) + ' doesn\'t contain reverse restriction site (BSMBI) \n'
    last = fragments[len(fragments)-1]
    if last.find(BSMBI_reverse) != (len(last)-len(BSMBI)):
        flags = flags + \
            'Last fragment doesn\'t contain restriction site (BSMBI) at its end \n'
    return flags

def filter_restriction_sgrna (df) :
    for index, row in df.iterrows() :
        test_restriction_sgrna=args.fragments[0] + row.Protospacer +args.fragments[1]
        if (test_restriction_sgrna.count(args.BSMBI) +  test_restriction_sgrna.count(BSMBI_reverse))>3:
            df=df.drop(index=index)

if __name__ == '__main__':
    args = argparser.parse_args()
    DataDir = os.path.abspath(os.path.dirname(__file__))
    flags = Assert_Fragments_BSMBI(args.BSMBI, args.fragments)
    if flags:
        raise Exception(flags)
    BSMBI_reverse = str(Seq(args.BSMBI).reverse_complement())
    primers = pd.read_csv(
        DataDir + '/data/sublibrary-primer-database-200.csv', index_col='Primer.Pair.Number')
    primers_reverse = primers.loc[args.Primer_pair]['Reverse.Binding']
    primers_forward = primers.loc[args.Primer_pair]['SKPP.forward.Seq']
    Study=pd.read_csv(args.Guide, sep=',',comment='#')
    Study.drop_duplicates(subset='Protospacer', inplace=True, keep=False)
    filter_restriction_sgrna(Study)
    guides=Study[['ID','Protospacer']]
    # Positive
    positive_options_df=pd.DataFrame([args.Vep[(0+i):(4+i)] for i in range(0,len(args.Vep),4)],columns = ['file','vep', 'N','Positive_Type'])
    positive_unused_df=pd.DataFrame([],columns=['ID','Protospacer'])
    for index, row in positive_options_df.iterrows():
        vep=pd.read_csv(row.vep, sep='\t',comment='#')
        Control_positive=vep[vep.Consequence==row.Positive_Type].Uploaded_variation
        csv=pd.read_csv(row.file,comment='#')
        Control_positive_df=csv.loc[[rowed.ID in list(Control_positive) for ind, rowed in csv.iterrows()]]
        Control_positive_df=Control_positive_df.loc[[rowed.ID not in guides.ID for ind, rowed in Control_positive_df.iterrows()]]
        filter_restriction_sgrna(Control_positive_df)
        positive_ran_df=Control_positive_df.sample(n=int(row.N),random_state=11)
        guides=pd.concat([guides,positive_ran_df[['ID','Protospacer']]])
        positive_unused_df=pd.concat([positive_unused_df,csv.loc[[w not in positive_ran_df.ID for w in csv.ID]]])
    #Negative
    negative_df=pd.DataFrame([args.Negative[(0+i):(2+i)] for i in range(0,len(args.Negative),2)],columns = ['file', 'N'])
    negative_unused_df=pd.DataFrame([],columns=['ID','Protospacer'])
    for index, row in negative_df.iterrows():
        csv = pd.read_csv(row.file, sep=',',comment='#')
        filter_restriction_sgrna(csv)
        csv=csv.loc[[rowed.ID not in guides.ID for ind, rowed in csv.iterrows()]]
        negative_ran_df=csv.sample(n=int(row.N),random_state=11)
        negative_unused_df=pd.concat([negative_unused_df,csv.loc[[w not in negative_ran_df.ID for w in csv.ID]]])
        guides=pd.concat([guides,negative_ran_df[['ID','Protospacer']]])
    guides=guides.drop_duplicates(subset='Protospacer', keep=False)
    SguidePerConcat = len(args.fragments)-1
    if len(guides.ID) % SguidePerConcat > 0 :
        df_unused=guides,negative_unused_df.loc[[w not in guides.Protospacer for w in negative_unused_df.Protospacer]]
        guides=pd.concat([guides,negative_unused_df.sample(n=SguidePerConcat-(len(guides.ID) % SguidePerConcat),random_state=11)])
    guides= guides.sample(frac=1,random_state=11)
    with open(args.Output, 'w') as out:
        for i in range(0, len(guides.index), SguidePerConcat):
            GuidesNames = []
            concat = [primers_forward]
            for j in range(0, SguidePerConcat):
                concat.extend([args.fragments[j], guides.iloc[i+j].Protospacer])
                GuidesNames.append(guides.iloc[i+j].ID)
            Names = ','.join(GuidesNames)
            concat.extend([args.fragments[SguidePerConcat],primers_reverse])
            concatemere = ''.join(concat)
            out.write('\t'.join([concatemere, Names]) + ' \n')
