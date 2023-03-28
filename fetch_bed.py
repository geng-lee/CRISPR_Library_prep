#!/usr/bin/python3

''' Crisper Library Design 1st Step

Goal provide a bed file with the sequences of all genes of interest.


'''

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os
import argparse
import re
import gffutils
import pyranges as pr

argparser = argparse.ArgumentParser(
    description='This software produces a bed file corresponding to the regions of interest as defined by a protein (corresponding to a genome) and the desired Feature. This represent the first step (potentially optional) in producing a library design crispr Array')
argparser.add_argument('-P', '--ProtFile', metavar='file', dest='prot_file',
                       type=str, required=True, help='List of protein to be examined')
argparser.add_argument('-G', '--Genome', metavar='name', dest='Genome', type=str, required=True,
                       help='prefix of Genome (*.fa)  and associated file (*.gff) to be used as reference (if custom feature used .gff not necessary)')
argparser.add_argument('-o', '--out', metavar='file', dest='Output',
                       type=str, required=False, default='out', help='Output')
argparser.add_argument('-B', '--Border', metavar='bp', dest='border', type=int, required=False,
                       default=30, help='Border to be included around feature as define in feature options')
argparser.add_argument('-F', '--Features', metavar='Annotation', dest='features', type=str,
                       required=False, nargs='*', default=['CDS'], help='Feature Types to be analysed (default exons)')
argparser.add_argument('-T', '--Tag', metavar='Annotation', dest='tag', type=str,
                       required=False, default='', help='Optional tag (usefull for Essential gene Controls)')
argparser.add_argument('-M','--Mane_Select', metavar = 'bool', dest = 'MANE_select_option', type = bool, required = False,default='True', help = 'Only select transcript with the main isoform (as define by MANE Select project)')



def Feature_Annotation(Feature, listed):
    returned = []
    if isinstance(Feature, str):
        h = list(
            set([y for x in q for y in db.children(x, featuretype=Feature)]))
    FeatureType = Feature
    for rec in h:
        Start_edit_freature = rec.start-args.border
        End_edit_freature = rec.end+args.border
        returned.append([rec.seqid, Start_edit_freature,
                         End_edit_freature])
    return returned


if __name__ == '__main__':
    args = argparser.parse_args()
    DataDir = os.path.abspath(os.path.dirname(__file__))
    prot = open(args.prot_file, 'r')
    Lines = prot.readlines()
    db = gffutils.FeatureDB("./data/genomes/"+args.Genome+'.gff.db')
    out = args.Output+'_'+args.tag+".bed"
    with open(out, "w") as output_handle:
        for line in Lines:
            protein = str(line.strip())
            if protein.startswith('custom:'):
                Custom = protein.split(':')
                Custom_rec = '\t'.join(
                    Custom[1], Custom[2], Custom[3], protein)
                output_handle.write()
            else:
                PAM_occurences = []
                if args.MANE_select_option:
                    q=[f.attributes['ID'][0] for f in db.children(db[protein].attributes['ID'][0]) if 'tag' in f.attributes.keys() and 'MANE_Select' in f.attributes['tag']]
                else:
                    q = [f.attributes['ID'][0]
                     for f in db.children(db[protein].attributes['ID'][0])]
                for Feature in args.features:
                    PAM_occurences.extend(Feature_Annotation(Feature, q))
                pyranges = pr.PyRanges(pd.DataFrame(PAM_occurences, columns=[
                                       'Chromosome', 'Start', 'End'])).merge()
                df = pyranges.as_df()
                df['protein'] = protein
                df.to_csv(output_handle, header=False, index=None, sep='\t')