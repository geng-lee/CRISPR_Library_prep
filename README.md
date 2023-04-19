# CRISPR_Library_prep
Script(s) and associated files for CRISPR guide RNA generation

Software necessary :

nextflow (https://www.nextflow.io/)
singulairty (https://docs.sylabs.io/guides/2.6/user-guide/quick_start.html)
python/3.8

Necessary Libraries :

(detailed in set-up.py)

Quick Start :

<ins> step 0) Software :  </ins>

Download this repository : git clone https://github.com/CERC-Genomic-Medicine/CRISPR_Library_prep.git (when public) 
Download CRISPOR : git clone https://github.com/maximilianh/crisporWebsite.git 
Install dependency : python3 set-up.py (preferably in a virtual environment see below) 

<ins> Step 1 Fetch_bed.py </ins> 

Goal : Create a merged region bed file corresponding to all specified protein feature 

Usage : python3 Fetch_bed.py -P your_list_protein.txt -G the_genome (additionnal options are availlable)  
 
| Option     | significance | Detail |
| ---      | ---       | ---       |
| -P | Protein list | File containing protein name (i.e. BRAC1) (1 per row) |
| -G     | Genome genome name | prefix of the genome (within directory data/genomes/[genome].gff3 |
| -o | out prefix | output prefix |
| -B | Border | Border to be included around feature as define in feature options (optional default '30')|
| -F | Feature type | Can be multiple (i.e -F 'CDS' 'Exon') (optional default 'CDS')|
| -M | gff record tagged by Mane (main isoforme project) | Flag (optional default false) |


<ins> Step 2 crispor.py </ins> 

Goal : Find and evaluate all guideRNA

Usage : python3 [PATH]/crispor.py genome --genomeDir=[PATH] [bedfile step1] [output.score] -o [output.offTarget] --maxOcc=6000 -p [PAM]

Detail see https://github.com/maximilianh/crisporWebsite/ for details

<ins> Step 3 Filter_offtarget.py </ins> 

Goal : This scrip allow for the filtering based on off-target activity prediction. 
 
python3 Filter_offtarget.py -I [output.offTarget] (additionnal options are availlable)   
 
argparser.add_argument('-I','--OffTarget', metavar = 'file', dest = 'OffTarget', type = str, required = True, help = 'Off Target file')
argparser.add_argument('-C','--CountTreshold', metavar = 'integer', dest = 'N', type = int, default=5, required = False, help = 'Off Target count threshold')
argparser.add_argument('-T','--Treshold', metavar = 'integer', dest = 'Thres', type = int, default=1.0, required = False, help = 'Threshold CFD to consider')
argparser.add_argument('-O','--out', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'Ouput file name')

| Option     | significance | Detail |
| ---      | ---       | ---       |
| -I | OffTarget file | OffTarget file produced by Crispor |
| -C     | CountTreshold | Off Target count threshold |
| -T | Treshold | Threshold CFD to consider (optional default: 1 |
| -O | Output | Ouput file name (optional default: 'out'|

 
<ins> Step 3 basic_annotation.py </ins> 

Goal : Create a file integrating editor information (possibly per-variant and per-Guide

python3 basic_annotation.py -E [editors] -S [output.score] -b [.bed step1] (additionnal options are availlable)  

| Option     | significance | Detail |
| ---      | ---       | ---       |
| -E | Editor | Editor type (will be examined against reference) (can be multiple i.e -E FLNS ABE8 |
| -S     | ScoreGuide | Crispor Score Guide Ouput |
| -O | out prefix | output prefix |
| -e | Exclude | list of guide to exclude (see step 3)|
| -V | per-Variant-VCF | flag (i.e usage : -V) to produce a per variant VCF suitable for VEP (optional ; default is False|
| -G | Per-Guide-VCF | flag (i.e usage : -G) to produce a per guide vcf file suitable for VEP (optional default False) |
| -L | length potospacer | length of the GuideRNA without PAM (optional ; default 20) |
| -B | bedfile | bedFile protein per region (see step 1) |

<ins> Step 4 VEP Annotation </ins> 

Goal : Annotate both the per-Variant-VCF and Per-Guide-VCF

See https://github.com/CERC-Genomic-Medicine/vep_pipeline for further information

usage : nextflow run util/VEP_variant/Annotation.nf -c nextflow.config   
  
<ins> Step 5 oligomer.py </ins>  
Goal : This script helps design oligomer/concatemere for array Crispr experiment
Usage: python3 oligomer.py -S [study.csv output step 3] -I [positive_control.csv output step 3] -V [positive_control.vep output step 4] -K [N istop] -N [negative_control.csv output step 3] -J [N negative_controls]  

| Option     | significance | Detail |
| ---      | ---       | ---       |
|  -S    |  Study File  | Study proteins/regions' csv (output of step 3) |
|   -I    |     IStop controls   | Positive Control proteins/regions' csv (output of step 3) |
|   -V    |     IStop controls VEP   | VEP of positive controls (output of step 4) |
|   -K    |     N Istop controls  | Number of Istop |
|   -N    |     Negative Controls   | Negative Control proteins/region's csv (output of step 3) |
|   -J    |     N Negative Controls   | Number of Negative Control associated Guides|
| -E | Exclude | list of guide to exclude (optional)|
| -O | out prefix | output prefix |
| -BSMBI | restriction site | Sequence of the restriction site (default 'CGTCTC')  |
| -P | Primer File | Primer file to be found in the data folder (example provided)|
| -f | Forward Primer index | Designation of the Primer pairs to be used |
| -F | Fragments | Oligomer fragment composition (optional ; default provided) (usage -F 'CGTCTCACACCG', 'GTTTTGAGACGgactgcCGTCTCcCACCG', 'GTTTaGAGACGggactaCGTCTCgCACCG', 'GTTTcGAGACGcttctcCGTCTCtCACCG', 'GTTTgGAGACG') |


