#! /bin/bash -l
#SBATCH -D /home/smhigdon/Projects/lactococcus 
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/ipscan-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/ipscan-stderr-%j.txt
#SBATCH -J ipscan
#SBATCH --array 5-6 # number of files to be processed
#SBATCH -t 120:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 16000
#SBATCH -p high
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu


# Name: ipscan_lactococcus.sh
# Created by: Shawn Higdon
# creation date: April 2, 2020
# This script will run InterproScan on nucleic acid fasta files generated from Prokka annotation of short read assemblies of individual microbial isolates.


# Load Modules

module load java/1.8
module load interproscan
source activate py27

OUTPUT_FOLDER=analysis

SEEDFILE=input_files/ips_inputs.txt

SEED=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

bname=`basename $SEED`
bname=`echo $bname | cut -d. -f1`

IPS_FOLDER=$OUTPUT_FOLDER/interproscan-bcw-lactococcus-only
mkdir -p $IPS_FOLDER

SAMPLE_FOLDER=$IPS_FOLDER/ips-${bname}
mkdir -p $SAMPLE_FOLDER

# Scan the nucleic acid fasta file using Interproscan v. 5.32-71.0

interproscan.sh -dp -t n -goterms -iprlookup -d $SAMPLE_FOLDER -f TSV,GFF3,HTML,SVG -i $SEED

source deactivate

hostname
