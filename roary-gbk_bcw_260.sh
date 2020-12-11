#!/bin/bash -l
#SBATCH -D /home/smhigdon/Projects/lactococcus
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/roary-all-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/roary-all-stderr-%j.txt
#SBATCH -J roary_all
#SBATCH -p high
#SBATCH -N 1	
#SBATCH -n 24
#SBATCH --mem=64000
#SBATCH -t 72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu

# Name: roary-bcw_lactococcus.sh
# Created by: Shawn Higdon
# creation date: June 30, 2020
# This script will generate a pan-genome for all lactococcus genomes using roary: the pan genome pipeline (v3.12.0 built with conda).

# Load Modules
source activate roary

# Define variables
ANALYSIS_FOLDER=analysis
mkdir -p $ANALYSIS_FOLDER

PROKKA_FOLDER=$ANALYSIS_FOLDER/prokka-gff-ms3

ROARY_FOLDER=$ANALYSIS_FOLDER/roary-gff-ms3
mkdir -p $ROARY_FOLDER

# Generate Pan-Genome with Roary using multi-fasta amino acid files from Prokka
roary \
	-e \
	-p 24 \
	-i 90 \
	-f $ROARY_FOLDER \
	-o all_isolate_clustered_proteins \
	$PROKKA_FOLDER/*.gff
	
source deactivate

hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks
stream
