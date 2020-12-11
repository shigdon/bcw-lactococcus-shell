#!/bin/bash -l
#SBATCH -D /home/smhigdon/Projects/lactococcus
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/scoary-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/scoary-rkr-stderr-%j.txt
#SBATCH -J scoary
#SBATCH -p high
#SBATCH -N 1	
#SBATCH -c 8
#SBATCH --mem=32000
#SBATCH -t 72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu

# Name: scoary_lactococcus.sh
# Created by: Shawn Higdon
# creation date: Mar 1, 2020
# This script will run Scoary: Microbial Pan-GWAS analysis on a pan-genome generated with roary: the pan genome pipeline (v3.12.0 built with conda).

# Load Modules
source activate roary

# Define variables
ANALYSIS_FOLDER=analysis
mkdir -p $ANALYSIS_FOLDER

ROARY_FOLDER=$ANALYSIS_FOLDER/roary-bcw-lactococcus-only
mkdir -p $ROARY_FOLDER

SCOARY_FOLDER=$ANALYSIS_FOLDER/scoary-bcw-lactococcus-only
mkdir -p $SCOARY_FOLDER

# Run Scoary on Pan-Genome from Roary 
scoary \
    -t input_files/scoary_bcw-lac-only_bnf.csv \
	-g $ROARY_FOLDER/gene_presence_absence.csv \
    -o $SCOARY_FOLDER \
    -p 1E-6 \
    -c BH \
    --upgma_tree \
    --threads 8 \
    --delimiter ,
	
source deactivate

hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks
stream
