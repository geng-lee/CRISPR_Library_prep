# CRISPR_Library_prep
Script(s) and associated files for CRISPR guide RNA generation

General steps :

1) Fetch_bed.py : Create a non-overlapping bed formated file of regions associated with specified proteins.
2) CRISPOR see below for installation.
2.5) Optional Filters
3) basic_annotation.py : Creates 1) reduced csv file of CRISPOR output linked to bedfile. 2) (potentially optional) vcf for VEP annotation 3)(optional) fetch all clinvar annotations associated with the CRISPOR guide RNA (for a specified editor(s))
4) (optional / necessary) VEP annotation.
5) Create concatemer libraries.

Software necessary :

nextflow singulairty
python2 python3

Necessary Libraries :

Bio  
pandas  
numpy  
argparse  
re  
gffutils  
pyranges as pr  

Quick Start :

Step 1) python3 Fetch_bed.py -P your_list_protein.txt -G the_genome (additionnal options are availlable)  
Step 2) python2 [PATH]/crispor.py genome --genomeDir=[PATH] [bedfile step1] [output.score] -o [output.offTarget] --maxOcc=6000 -p [PAM]  
Step 2.5) python3 Filter_offtarget.py -I [output.offTarget] (additionnal options are availlable)  
Step 3) python3 basic_annotation.py -E [editors] -S [output.score] -b [.bed step1] (additionnal options are availlable)  
step 4) nextflow run util/VEP_variant/Annotation.nf -c nextflow.config   
  
... To produce a full Crispr library experiment positive controls (induced stop) and negative controls are necessary.  
  
step 5) python3 oligomer.py -S [study.csv output step 3] -I [positive_control.csv output step 3] -V [positive_control.vep output step 4] -K [N istop] -N [negative_control.csv output step 3] -J [N negative_controls]  


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
