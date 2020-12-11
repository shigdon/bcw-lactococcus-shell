#!/bin/bash -l
#SBATCH -D /home/smhigdon/Projects/lactococcus
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/prok_bcw-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/prok_bcw-stderr-%j.txt
#SBATCH -J prokka_bcw
#SBATCH --array 12-14#number of files to be processed
#SBATCH -p high
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=8000
#SBATCH -t 72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu

# Name: prokka_all_mucilage_genomes_array.sh
# Created by: Shawn Higdon
# creation date: Mar 1, 2020
# This script will annotate short read assemblies from individual bacterial genomes using the PROKKA pipeline for BCW mucilage isolate genomes - as an array job.

# Load Modules
source activate prokka

# Define variables
OUTPUT_FOLDER=analysis
mkdir -p $OUTPUT_FOLDER

SEEDFILE=input_files/prokka_inputs.txt
# bcw_prokka_inputs.txt is the list of files you want to process, 1 per line

SEED=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
bname=`basename $SEED`
sample=$(echo $bname | cut -d. -f1)

SAMPLE_FOLDER=$OUTPUT_FOLDER/$sample
mkdir -p $SAMPLE_FOLDER

# Annotate the MEGAHIT final.contigs.fa file using PROKKA
PROKKA_FOLDER=$SAMPLE_FOLDER/prokka

prokka \
	--outdir $PROKKA_FOLDER \
	--prefix $sample \
	--locustag $sample \
	--cpus 8 \
    --force \
	$SEED

source deactivate

hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks
stream
