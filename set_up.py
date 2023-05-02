#!/usr/bin/python3

## Use python/3.8


import sys
import subprocess
import pkg_resources
import os

required = {'wget', 'gffutils','Bio', 'pandas','pytabix','matplotlib','numpy','argparse','pyranges','biopython','numpy','scikit-learn','twobitreader','xlwt','rs3'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    python = sys.executable
    subprocess.check_call(['pip3', 'install', *missing], stdout=subprocess.DEVNULL)

subprocess.check_call(['git', 'clone', 'https://github.com/maximilianh/crisporWebsite.git'], stdout=subprocess.DEVNULL)

#### Download Accessory files
import wget
import gffutils


DataDir=os.path.abspath(os.path.dirname(__file__))
os.chdir(DataDir+'/data/genomes')

url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.basic.annotation.gff3.gz'
gff3 = wget.download(url)
os.rename('gencode.v42.basic.annotation.gff3.gz', 'hg38.gff3.gz')

def transform_func(x):
    # adds some text to the end of transcript IDs
    if 'transcript_id' in x.attributes:
        x.attributes['transcript_id'][0] += '_transcript'
    return x

# db = gffutils.create_db('hg38.gff3.gz', id_spec={'gene': 'gene_id', 'transcript': "transcript_id"}, merge_strategy="create_unique", transform=transform_func, keep_order=True)


url='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz'
wget.download(url)

url='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi'
wget.download(url)

#cd ../../
#mkdir genomes
#cd genomes
#mkdir hg38
#cd hg38
#wget -r -l1 --no-parent -nd  --reject 'index*' --reject 'robots*' http://crispor.tefor.net/genomes/hg38/
