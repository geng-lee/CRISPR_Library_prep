# CRISPR_Library_prep
Script(s) and associated files for CRISPR guide RNA generation

General steps :

1) Fetch_bed.py : Create a non-overlapping bed formated file of regions associated with specified proteins.
2) CRISPOR see below for installation.
3) basic_annotation.py : Creates 1) reduced csv file of CRISPOR output linked to bedfile. 2) (potentially optional) vcf for VEP annotation 3)(optional) fetch all clinvar annotations associated with the CRISPOR guide RNA (for a specified editor(s))
4) (optional / necessary) VEP annotation.
5) Create concatemer libraries.

Necessary Libraries :

Bio  
pandas  
numpy  
argparse  
re  
gffutils  
pyranges as pr  


Crispor installation on ComputeCanada :

module load python/2.7.18  
python -m venv $HOME/env_CRISPOR  
source $HOME/env_CRISPOR/bin/activate  
module load mysql  
module load blat/3.5  
module load bwa kentutils/401 gffread/0.12.3  
pip install scikit-learn==0.16.1 pandas twobitreader biopython xlwt  
pip install matplotlib  
pip install MySQL-python  
https://github.com/maximilianh/crisporWebsite.git  
