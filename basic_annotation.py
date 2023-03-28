#!/usr/bin/python3

'''
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 1.0 (beta)
YEAR: 2023

Description Pipeline

1) Get BedFile from Protein List (Optional)
2) Run CRISPOR
3) Run CDF_Filter_script (Optional)
4) Run Basic_Annotation_script.py (Clinvar (Optional) / Pre-VEP vcf File (Optional) / General)
5) Run VEP (Optional / Necessary for istop controls)
6) Run OLIGOMER_prep.py

'''


from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os
import argparse
import re
import gffutils
import pysam
import gzip
import pyranges as pr

argparser = argparse.ArgumentParser(description = 'This Software is a part of a pipeline tool to design Guiding RNA. \n This script relates back the CRISPOR output an initial genome annotation and optionnaly can provide editor specific annotation (Clinvar and/or VCF file)')
argparser.add_argument('-E','--Editor', metavar = 'file name', dest = 'Editor', nargs='+', type = str, required = False, default=None, help = 'Editor type (will be examined against reference)')
argparser.add_argument('-S','--ScoreGuide', metavar = 'file name', dest = 'scoreGuide', type = str, required = True, help = 'Crispor score ouput')
argparser.add_argument('-O','--Output', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'prefix of outputs')
argparser.add_argument('-e','--Exclude', metavar = 'file', dest = 'excludes', type = str, required = False, default ='', help = 'List of Guides to exclude')
argparser.add_argument('-C','--Clinvar', dest = 'clinvar', action='store_true', help = 'flag to analyse against clinvar')
argparser.add_argument('-V','--VEP', dest = 'VEP', action='store_true', help = 'flag to produce vcf file suitable for VEP')
argparser.add_argument('-L','--length', metavar = 'int', dest = 'length', type = int, required = False, default ='20', help = 'length of the GuideRNA without PAM')
argparser.add_argument('-B','--bed', metavar = 'file', dest = 'bed', type = str, required = True, help = 'bedFile protein per region')



if __name__ == '__main__':
        args = argparser.parse_args()
        DataDir=os.path.abspath(os.path.dirname(__file__))
        if not os.path.exists(DataDir+'/data/Editor_library.csv') :
                raise Exception('Necessary auxiliairy files are absent.')
        else :
                if not args.Editor == None :
                        data=pd.read_csv(DataDir+'/data/Editor_library.csv',index_col='Editor',sep='\t')
                        editor=data.loc[args.Editor]
                        del data
        scoreGuides=pd.read_csv(args.scoreGuide,sep='\t')
        scoreGuides.drop_duplicates(subset='targetSeq', inplace=True, keep=False)
        scoreGuides['targetSeq_plusStrand']=[str(Seq(row.targetSeq).reverse_complement()) if 'rev' in  row['guideId']  else row['targetSeq'] for index, row in scoreGuides.iterrows()]
        scoreGuides['protospacer']=[j[0:args.length].upper() for j in scoreGuides.targetSeq]
        scoreGuides['PAM']=[j[args.length:len(j)].upper() for j in scoreGuides.targetSeq]
        scoreGuides.drop_duplicates(subset='targetSeq_plusStrand', inplace=True, keep=False)
        if not args.excludes == "" :
                excludes=pd.read_csv(args.excludes,sep='\t')
                scoreGuides=scoreGuides.loc[scoreGuides['targetSeq'].isin(excludes[0].tolist())]
        ### Find protein corresponding
        scoreGuides['positions']=[int(i.replace('rev','').replace('forw','')) for i in scoreGuides.guideId]
        scoreGuides['Chromosome']=[re.split('-|:', j)[0] for j in scoreGuides['#seqId']]
        scoreGuides['start_seqId']=[int(re.split('-|:', j)[1]) for j in scoreGuides['#seqId']]
        scoreGuides['target_Start']=scoreGuides['start_seqId']+scoreGuides['positions']
        ### Find proiten
        bed=pr.read_bed(args.bed)
        protein=[]
        for index, row in scoreGuides.iterrows():
                df=pd.DataFrame({"Chromosome": [row['Chromosome']], "Start": [row.target_Start-args.length if 'rev' in  row.guideId  else row.target_Start ],"End":[row.target_Start if 'rev' in  row.guideId  else row.target_Start + args.length]})
                position= pr.PyRanges(df)
                names=[i for i in position.join(bed).Name]
                if len(set(names))>1:
                        raise Exception('One protein overlaps with anoter.')
                else :
                        protein.extend(list(set(names)))
        scoreGuides['Protein']=protein
        scoreGuides['strand']=['-'if 'rev' in j  else '+' for j in scoreGuides.guideId]
        scoreGuides.sort_values(['Protein','guideId','strand'])
        scoreGuides['num']= scoreGuides.groupby(['Protein']).cumcount()+1
        scoreGuides['ID']=scoreGuides['Protein']+'_'+ [str(j) for j in scoreGuides['num']]
        scoreGuides.index=scoreGuides['ID']
        with  open(args.Output+'_general.csv', 'wt') as general :
                general.write('ID,gRNA_seq,gRNA_seq_POSstrand,Chromosome,POSstart,POSend,strand,cfdSpecScore,Doench\'16-Score\n')
                for index, row in scoreGuides.iterrows():
                        general.write('{ID},{protospacer},{gRNA_seq},{gRNA_seq_POS},{chrom},{POSstart},{Strand},{cfdSpecScore},{Doench}\n'.format(
                                ID=row.ID,
                                protospacer=row['protospacer'],
                                gRNA_seq=row.targetSeq,
                                gRNA_seq_POS=str(Seq(row.targetSeq).reverse_complement()) if 'rev' in  row['guideId']  else row['targetSeq'],
                                chrom=row['Chromosome'],
                                Strand=row.strand,
                                POSstart=str(row.target_Start),
                                cfdSpecScore=str(row.cfdSpecScore),
                                Doench=str(row['Doench \'16-Score'])))
        ### producting Edditor specific files
        if not args.Editor == None :
                for index, i in editor.iterrows() : ###***?!!?***
                        editor_df=scoreGuides.copy()
                        i.window_start=i.window_start
                        i.window_end=i.window_end
                        editor_df['editing_windowSeq']=[j[i.window_start-1:i.window_end].upper() for j in editor_df.targetSeq]
                        editor_df['editing_windowSeq']= [str(Seq(editor_df['editing_windowSeq'][j]).reverse_complement()) if 'rev' in  editor_df['guideId'][j]  else editor_df['editing_windowSeq'][j]   for j in editor_df.index]
                        editor_df['is_empty'] = [not (editor.BE[0] in x.upper() or editor.BE_RC[0] in x.upper()) for x in editor_df["editing_windowSeq"]]
                        editor_df['editing_windowSTART']=[editor_df['target_Start'][j]-args.length+i.window_start-1 if 'forw' in editor_df.guideId[j] else editor_df['target_Start'][j] +len(editor_df.targetSeq[j])-i.window_end for j in editor_df.index]
                        editor_df['editing_windowEND']=editor_df['editing_windowSTART'] + (i.window_end-i.window_start)
                        editor_df['editing_window_mutated']=[row.editing_windowSeq.replace(i.BE[0],i.BE[1]) if row.strand=='+' else row.editing_windowSeq.replace(i.BE_RC[0],i.BE_RC[1]) for indexed,row in editor_df.iterrows()]
                        editor_df['nchange']=[sum([not rowed.editing_windowSeq[j]==rowed.editing_window_mutated[j] for j in range(0,len(rowed.editing_window_mutated))]) for indexed, rowed in editor_df.iterrows()]
                        with open(args.Output+'.'+i.name +'_empties' + '.txt', 'w') as empties:
                                for index, row in editor_df.iterrows():
                                        if row.editing_windowSeq ==row.editing_window_mutated :
                                                empties.write(str(row.name) +'\n')
                        if args.VEP :
                                with gzip.open(args.Output+'.'+i.name+'.vcf.gz','wt') as vcf :
                                        vcf.write('##fileformat=VCFv4.1\n')
                                        vcf.write('##CrisprTransition='+str(index)+ '\n')
                                        vcf.write('##contig=<ID=chr1,length=249250621>\n##contig=<ID=chr10,length=135534747>\n##contig=<ID=chr11,length=135006516>\n##contig=<ID=chr12,length=133851895>\n##contig=<ID=chr13,length=115169878>\n##contig=<ID=chr14,length=107349540>\n##contig=<ID=chr15,length=102531392>\n##contig=<ID=chr16,length=90354753>\n##contig=<ID=chr17,length=81195210>\n##contig=<ID=chr18,length=78077248>\n##contig=<ID=chr19,length=59128983>\n##contig=<ID=chr2,length=243199373>\n##contig=<ID=chr20,length=63025520>\n##contig=<ID=chr21,length=48129895>\n##contig=<ID=chr22,length=51304566>\n##contig=<ID=chr3,length=198022430>\n##contig=<ID=chr4,length=191154276>\n##contig=<ID=chr5,length=180915260>\n##contig=<ID=chr6,length=171115067>\n##contig=<ID=chr7,length=159138663>\n##contig=<ID=chr8,length=146364022>\n##contig=<ID=chr9,length=141213431>\n##contig=<ID=chrM,length=16571>\n##contig=<ID=chrX,length=155270560>\n##contig=<ID=chrY,length=59373566> \n')
                                        vcf.write('##INFO=<ID=Protospacer,Number=.,Type=String,Description="Protospace"> \n')
                                        vcf.write('##INFO=<ID=STRAND,Number=.,Type=String,Description="Strand"> \n')
                                        vcf.write('##INFO=<ID=PAM,Number=.,Type=String,Description="Protospacer adjacent motif"> \n')
                                        vcf.write('##INFO=<ID=Nchange,Number=.,Type=Integer,Description="Nummber of nucleotide changes"> \n')
                                        vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                                        for index, row in editor_df.iterrows():
                                                if not row.editing_windowSeq ==row.editing_window_mutated :
                                                        vcf.write('{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{END}\t{FILTER}\t{INFO}\n'.format(
                                                                chrom = str(row['Chromosome']),
                                                                pos = str(row['editing_windowSTART']) ,
                                                                vid=str(row.ID),
                                                                ref= str(row['editing_windowSeq']) ,
                                                                alt= str(row['editing_window_mutated']) ,
                                                                END='.',
                                                                FILTER ='.' ,
                                                                INFO='Protospacer='+row['protospacer']+ ';PAM='+ row['PAM'] +';STRAND='+ row.strand + ';Nchange=' + str(row.nchange)))
                        if args.clinvar:
                                with open(args.Output+'.'+i.name+'.Clinvar.csv', "w") as out :
                                        vcf=pysam.VariantFile(DataDir+'/data/genomes/clinvar.vcf.gz')
                                        clinvar_version = None
                                        for line in vcf.header.records:
                                                if line.key == 'fileDate':
                                                        clinvar_version = line.value
                                                        break
                                        out.write('## Clinvar file Date : ' + clinvar_version + '\n')
                                        out.write(",".join(['Chromosome','Position','Reference','Alternate','GuideRNA_ID','GuideRNA_sequence','CLNSIG','CLNDN','MC'])+'\n')
                                        ranges=editor_df[['Chromosome','editing_windowSTART','editing_windowEND','strand']]
                                        ranges.rename(columns={'Chromosome':"Chromosome","editing_windowSTART":'Start', 'editing_windowEND':'End'}, inplace=True)
                                        ranges['Start']=ranges.Start
                                        ranges['End']= ranges.End +1      # To include the ending
                                        ranges.index=range(0,len(ranges))
                                        ranges['Sequence']=[j for j in editor_df.protospacer]
                                        ranges['strand']=[j for j in editor_df.strand]
                                        ranges['ID']=[j for j in editor_df.ID]
                                        unmerged=pr.PyRanges(ranges)
                                        merged=unmerged.merge().as_df()
                                        merged.to_csv(args.Output+'.'+i.name+'.bed', index=False)
                                        for rangeIndex, RangeRow in merged.iterrows() :
                                                for rec in vcf.fetch(RangeRow['Chromosome'].strip('chr'), RangeRow.Start, RangeRow.End-1) : ## pysam already includes ending
                                                        if not rec.alts is None:
                                                                for alt in rec.alts :
                                                                        if len(rec.ref)==1 and (rec.ref==i.BE[0] and i.BE[1] == alt) :
                                                                                df=pd.DataFrame({"Chromosome": ['chr'+str(rec.contig)], "Start": [rec.pos],"End":[rec.pos+1]})
                                                                                position= pr.PyRanges(df)
                                                                                overlap=position.join(unmerged).as_df()
                                                                                overlap=overlap[overlap.strand=='+']
                                                                                CLNDN=[""] if rec.info.get('CLNDN') is None else rec.info.get('CLNDN')
                                                                                CLNSIG=[""] if rec.info.get('CLNSIG') is None else rec.info.get('CLNSIG')
                                                                                if not overlap.empty:
                                                                                        out.write('chr'+str(rec.contig) + ','+str(rec.pos) + ',' + rec.ref + ','+ alt+ ','+'|'.join(overlap.ID)+',' + '|'.join(overlap.Sequence) + ',' + '|'.join([j for j in CLNSIG]) + ',' + '|'.join([j for j in CLNDN])+ ','+ "|".join([j for j in rec.info.get('MC')])+'\n' )
                                                                        if len(rec.ref)==1 and (rec.ref==i.BE_RC[0] and i.BE_RC[1] == alt):
                                                                                df=pd.DataFrame({"Chromosome": ['chr'+str(rec.contig)], "Start": [rec.pos],"End":[rec.pos+1]})
                                                                                position= pr.PyRanges(df)
                                                                                overlap=position.join(unmerged).as_df()
                                                                                overlap=overlap[overlap.strand=='-']
                                                                                CLNDN=[""] if rec.info.get('CLNDN') is None else rec.info.get('CLNDN')
                                                                                CLNSIG=[""] if rec.info.get('CLNSIG') is None else rec.info.get('CLNSIG')
                                                                                if not overlap.empty:
                                                                                        out.write('chr'+str(rec.contig) + ','+str(rec.pos) + ',' + rec.ref + ','+ alt+ ','+'|'.join(overlap.ID)+',' + '|'.join(overlap.Sequence) + ',' + '|'.join([j for j in CLNSIG]) + ',' + '|'.join([j for j in CLNDN])+ ','+ "|".join([j for j in rec.info.get('MC')])+'\n' )

                                vcf.close()
####################################################
#for index, row in clinvar.iterrows() :
#...     splited=row.GuideRNA.split('|')
#...     for j in splited :
#...             if guides.loc[j].POS>row.Position or guides.loc[j].POS+6<row.Position:
#...                     print(row)
#...                     print(j)
#...                     break