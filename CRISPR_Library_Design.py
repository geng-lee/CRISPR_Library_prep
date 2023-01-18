#!/usr/bin/python3

''' Crisper Library Design : 

This Software is a tool to design Guiding RNA.
This tool is flexible, permitting control design, clinvar annotation, 
annotated and custom feature.

Usage :

python3 CRISPR_Library_Design.py -P protein_file -O Ortholog -E Editor -G genome -F exon

Control options:

--empty : produces sgRNA without editor action
--Feature with options start_codon | stop_codon
--ProteinFile with Essential or non-essential Gene

'''

'''
TO DO : 
--Integration of CFD score calculation (based on Doench et al 2016)
--Integration of Oligo Design
--De novo Stop/Start codon analysis

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


argparser = argparse.ArgumentParser(description = 'This Software is a tool to design Guiding RNA. \n This tool is flexible, permitting control design, clinvar annotation, \n annotated and custom feature. ')
group = argparser.add_mutually_exclusive_group()
argparser.add_argument('-P','--ProtFile', metavar = 'file', dest = 'prot_file', type = str, required = True, help = 'List of protein to be examined')
argparser.add_argument('-O','--Ortholog', metavar = 'name', dest = 'Type', type = str, required = True, help = 'CRISPER type (will be examined against reference)')
argparser.add_argument('-E','--Editor', metavar = 'name', dest = 'Editor', type = str, required = True, help = 'Editor type (will be examined against reference)')
argparser.add_argument('-G','--Genome', metavar = 'name', dest = 'Genome', type = str, required = True, help = 'prefix of Genome (*.fa)  and associated file (*.gff3) to be used as reference (if custom feature used .gff3 not necessary)')
argparser.add_argument('-o','--out', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'Output')
argparser.add_argument('-B','--Border', metavar = 'name', dest = 'border', type = int, required = False, default=30, help = 'Bordering region to be included around feature as define in feature options (default 30)')
argparser.add_argument('-C','--clinvar', dest = 'clinvar', required = False, action='store_false', help = 'Whether or not to annotate feature with clinvar (default:True ; Requieres clinvar associated files)')
argparser.add_argument('-e','--empty', dest = 'empty', required = False, action='store_true', help = 'Produces only ')
argparser.add_argument('-l','--limit', dest = 'limit',type=int,required = False, default=0, help = 'Editor type (will be examined against reference)')
group.add_argument('-F','--Features', metavar = 'Annotation', dest = 'features',type=str, required = False, nargs='*',default=['exon'], help = 'Feature Types to be analysed (default exons)')
group.add_argument('-f','--Custom_feature', metavar = 'Region', dest = 'Custom_feature',type=str, required = False, nargs=1, help = 'Optional : Feed a custom Region (style chrom:beg:end)')

### Finder concept
# 1) Fech actual sequence
# 2) Look for PAM or PAM_reverse complement
# 3) Ends (+) or Starts (-) are essentially PAM location not relative (so ajusted for begining of exon)
# 4) define oposite of Starts/Ends
# 5) get sequence in between start and end on correct strand


def Finder(genome,chrom,beg,end,strand) :
	seq=str(genome[chrom].seq[beg:end])
	if strand=='+':
		starts=np.add(beg,[m.start() for m in re.finditer(crisper.PAM, str(seq))]).tolist()
		ends=np.add(crisper.gRNA_Length+len(crisper.PAM),starts).tolist()
		sequence=[str(genome[chrom].seq[starts[i]:ends[i]]) for i in range(0,len(starts))]
		editing_window=[i[len(crisper.PAM)+editor.window_start:len(crisper.PAM)+editor.window_end] for i in sequence]
		returned=[starts,ends,['+' for i in range(0,len(starts))],sequence,[crisper.PAM for i in range(0,len(starts))],editing_window]
	if strand=='-':
		ends=np.add(beg,[m.end() for m in re.finditer(crisper.PAM_RC, str(seq))]).tolist()
		starts=np.subtract(ends,crisper.gRNA_Length+len(crisper.PAM)).tolist()
		sequence=[str(genome[chrom].seq[(starts[i]):(ends[i])]) for i in range(0,len(starts))]
		editing_window=[i[len(i)-len(crisper.PAM)-editor.window_end:len(i)-len(crisper.PAM)-editor.window_start] for i in sequence]
		returned=[starts,ends,['-' for i in range(0,len(starts))],sequence,[crisper.PAM for i in range(0,len(starts))],editing_window]
	return returned


def Feature_Annotation(Feature) :
	returned=[]
	if isinstance(Feature,str):
		h=list(set([y for x in q for y in db.children(x,featuretype=Feature)]))
		FeatureType=Feature
	else :
		h=[Feature]
		FeatureType='Custom'
	for rec in h:
		for strand in ['+','-']:
			Start_edit_freature=(rec.start-len(crisper.PAM)-editor.window_end, rec.start+len(crisper.PAM)+editor.window_start)[strand=='-']
			End_edit_freature=(rec.end-len(crisper.PAM)-editor.window_start, rec.end+len(crisper.PAM)+editor.window_end)[strand=='-']
			Finds=Finder(genome,rec.seqid,Start_edit_freature-args.border,End_edit_freature+args.border,strand)
			occ=len(Finds[0])
			identifier=[[protein for x in range(0,occ)],[FeatureType for x in range(0,occ)],[rec.seqid for x in range(0,occ)],[rec.strand for x in range(0,occ)],[rec.start for x in range(0,occ)],[rec.end for x in range(0,occ)]]
			if Finds :
				if not returned :
					returned=identifier+Finds
					Finds=None
				else : 
					returned=[i + j for i,j in zip(returned, identifier+Finds)] 
					Finds=None
	if not h :
		raise Exception(protein + ' does not possess ' + Feature + ' in corresponding GTF file')
	return returned

if __name__ == '__main__':
  ### Loading and associated test
	args = argparser.parse_args()
	DataDir=os.path.abspath(os.path.dirname(__file__))
	if not os.path.exists(DataDir+'/data/crisper_orthologs.csv') :
		raise Exception('Necessary auxiliairy files are absent.')
	else :
		data=pd.read_csv(DataDir+'/data/crisper_orthologs.csv',index_col='Cas_ortholog')
		crisper=data.loc[args.Type]
		del data
	if not os.path.exists(DataDir+'/data/Editor_library.csv') :
		raise Exception('Necessary auxiliairy files are absent.')
	else :
		data=pd.read_csv(DataDir+'/data/Editor_library.csv',index_col='Editor',sep='\t')
		editor=data.loc[args.Editor]
		del data
	if args.clinvar:
		if not os.path.exists(DataDir+'/data/genomes/Clinvar.vcf.gz') or not os.path.exists(DataDir+'/data/genomes/Clinvar.vcf.gz.tbi') :
			raise Exception('Clinvar File are absent.')
	if not args.Custom_feature:
		if not os.path.exists("./data/genomes/"+args.Genome+'.gff.db') :
				raise Exception('Annotation file is absent despite (only allowed if look-up is custom feature)')
		else:
			db = gffutils.FeatureDB("./data/genomes/"+args.Genome+'.gff.db')
### Look for Annotations and genome
	prot= open(args.prot_file, 'r')
	Lines = prot.readlines()
	db = gffutils.FeatureDB("./data/genomes/"+args.Genome+'.gff.db')
	genome = SeqIO.index("./data/genomes/"+args.Genome+'.fa','fasta')
	PAM_occurences = []
	### look for guide RNA
	if not args.Custom_feature:
		for line in Lines:			## For each protein
			protein=line.strip()
			q=[f.attributes['ID'][0] for f in db.children(db[protein].attributes['ID'][0])]
			for Feature in args.features:
				if not PAM_occurences:
					PAM_occurences=Feature_Annotation(Feature)
				else :
					PAM_occurences=[i + j for i,j in zip(PAM_occurences,Feature_Annotation(Feature))]
	if args.Custom_feature:
			print(args.Custom_feature[0])
			split_costum_feature=args.Custom_feature[0].split(':')
			protein='Custom_feature'
			Feature=gffutils.Feature(id='Custom_feature',seqid=split_costum_feature[0],start=split_costum_feature[1],end=split_costum_feature[2],strand='None')
			PAM_occurences=Feature_Annotation(Feature)
	DF=pd.DataFrame(PAM_occurences).T
	DF.columns=['protein','feature','chromosome','feature_strand','feature_start','feature_end','SG_start', 'SG_end','SG_strand','SG_sequence','PAM_sequence','editing_window']
	DF.drop_duplicates('SG_sequence',keep='first',inplace=True)
	print(DF)
	DF.sort_values(by=["protein", "SG_start"], inplace=True)
	with open(args.Output, 'w') as out:
		if args.empty:
			is_empty = [not (editor.BE[0] in x.upper() or editor.BE_RC[0] in x.upper()) for x in DF["editing_window"]]
			output=DF.iloc[is_empty]
			if args.limit==0:
				output.to_csv(out,sep='\t',index=False)
			else :
				output.sample(n=args.limit,axis=0).to_csv(out,sep='\t',index=False)
		else :
			is_not_empty = [(editor.BE[0] in x.upper() or editor.BE_RC[0] in x.upper()) for x in DF["editing_window"]]
			DF=DF.iloc[is_not_empty]
			if not args.limit==0:
				DF=DF.sample(n=args.limit,axis=0)
			if args.clinvar: 
				vcf=pysam.VariantFile(DataDir+'/data/genomes/Clinvar.vcf.gz')
			for rows in range(0,len(DF)-1):
				row=DF.iloc[rows]
				beg=(row['SG_start']+editor.window_start + len(crisper.PAM),row['SG_end']-editor.window_end - len(crisper.PAM))[row['SG_strand']=='-']
				for pos in range(0,(editor.window_end-editor.window_start)) :
					if (beg+pos)<=(row['feature_end']+args.border) and (beg+pos)>=(row['feature_start']+args.border):
						SGref=str(genome[row['chromosome']].seq[beg+pos-1]).upper()
						print(str(beg+pos))
						print(SGref)
						if SGref==editor.BE[0] or editor.BE_RC[0]==SGref:
							empty=True
							if args.clinvar :
								for rec in vcf.fetch(row['chromosome'].strip('chr'), beg+pos-1, beg+pos+1) :
									if rec.alts is not None: ### Extremely odd occurence 
										if rec.ref==SGref and (editor.BE[1] in rec.alts or editor.BE_RC[1] in rec.alts):
											empty=False
											out.write('\t'.join(str(i) for i in row) + '\t'+str(beg+pos) + '\t' + SGref + '_'+ (editor.BE_RC[1],editor.BE[1])[SGref==editor.BE[0]]+ '\t' + str(rec.info.get('CLNSIG')) + ' ' + str(rec.info.get('CLNDN'))+ ' '+ str(rec.info.get('MC'))+'\n' )
							if empty :
								out.write('\t'.join(str(i) for i in row) + '\t'+str(beg+pos) + '\t' + SGref+(editor.BE_RC[1],editor.BE[1])[SGref==editor.BE[0]] +'\t' + ' \n' )


