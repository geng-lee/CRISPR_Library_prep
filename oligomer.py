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
argparser.add_argument('-I', '--iStop', metavar='name',
                       dest='iStop', type=str, required=True, help='iStop Positive Controls')
argparser.add_argument('-V', '--iStopVep', metavar='name',
                       dest='Vep', type=str, required=True, help='Variant effect predictor associated with Positive control library')
argparser.add_argument('-K', '--N_iStop', metavar='integer',
                       dest='NiStop', type=int, required=True, help='Number of Positive controls')
argparser.add_argument('-N', '--Negative', metavar='name',
                       dest='Negative', type=str, required=True, help='Negative controls ')
argparser.add_argument('-J', '--N_Negative', metavar='integer',
                       dest='N_negative', type=int, required=True, help='number of negative controls')
argparser.add_argument('-E','--Exclusion', metavar = 'file', dest = 'Exclusion', type = str, required = False,default=None, help = 'List of guides to be excluded')
argparser.add_argument('-O', '--out', metavar='file', dest='Output',
                       type=str, required=False, default='out', help='Ouput file name')
argparser.add_argument('--BSMBI', metavar='file', dest='BSMBI',
                       type=str, required=False, default='CGTCTC', help='Restriction Site')
argparser.add_argument('-P', '--PRIMER_FILE', metavar='file',
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
        test_restriction_sgrna=args.fragment[0] + row.gRNA_seq +args.fragment[1]
            if (test_restriction_sgrna.count(args.BSMBI) +  test_restriction_sgrna.count(BSMBI_reverse))>3:
                df.drop(index=index, inplace=True)

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
    # process CRISPOR file
    vep=pd.read_csv(args.Vep, sep='\t',index_col='#Uploaded_variation',skiprows=90)
    Stoppers=vep[vep.Consequence=='stop_gained'].index
    istop = pd.read_csv(args.iStop, sep=',',index_col='ID').loc[Stoppers].sample(n=args.NiStop,random_state=11)
    negative = pd.read_csv(args.Negative, sep=',',index_col='ID').sample(n=args.N_negative,random_state=11)
    guides=pd.read_csv(args.Guide, sep=',',index_col='ID')
    guides.drop_duplicates(subset='gRNA_seq_POSstrand', inplace=True, keep=False)
    guides=pd.concat([guides,istop,negative])
    if args.Exclusion :
        exclusion=pd.read_csv(args.Exclusion, sep='\t')
        guides=guides.loc[guides['targetSeq'].isin(excludes[0].tolist())]
    filter_restriction_sgrna(guides)
    guides= guides.sample(frac=1,random_state=11)
    SguidePerConcat = len(args.fragments)-1
    with open(args.Output, 'w') as out:
        for i in range(0, len(guides.index) // SguidePerConcat):
            GuidesNames = []
            concat = [primers_forward]
            for j in range(0, SguidePerConcat):
                concat.extend([args.fragments[j], guides.iloc[i * SguidePerConcat+j].targetSeq])
                GuidesNames.append(guides.iloc[i*SguidePerConcat+j].label)
            Names = ','.join(GuidesNames)
            concat.extend([args.fragments[SguidePerConcat],primers_reverse])
            concatemere = ''.join(concat)
            out.write('\t'.join([concatemere, Names]) + ' \n')
