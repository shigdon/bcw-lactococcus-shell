#!/bin/bash -l
#SBATCH -D /home/smhigdon/Projects/lactococcus
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/bwa_lactococcus-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/bwa_lactococcus-stderr-%j.txt
#SBATCH -J bwa_lacto
#SBATCH --array 3-72 # number of files to be processed
#SBATCH -p high
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8000
#SBATCH -t 96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu

# Name: bwa_lactococcus.sh
# Created by: Shawn Higdon
# creation date: June 22, 2020
# Description: A script to map trimmed fastq paired end reads to genome assemblies from MEGAhit

# Load modules
module load bwa

# Define Variable

OUTPUT_FOLDER=analysis # folder for all output
mkdir -p $OUTPUT_FOLDER

SEEDFILE=input_files/bwa_asm_inputs.txt # the list of files that will be processed (absolute path), 1 per line.

SEED=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

bname=`basename $SEED`
bname=`echo $bname | cut -d. -f1`

SAMPLE_FOLDER=$OUTPUT_FOLDER/$bname
mkdir -p $SAMPLE_FOLDER

MAPPING_FOLDER=$SAMPLE_FOLDER/bwa_mapping
mkdir -p $MAPPING_FOLDER

# Step One: index the assemblies with bwa version 0.7.9a

bwa index $SEED

# Step Two: map the reads with bwa

bwa mem $SEED $SAMPLE_FOLDER/*.trim_R1.fq.gz $SAMPLE_FOLDER/*.trim_R2.fq.gz > $MAPPING_FOLDER/${bname}.aln.sam


hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks
stream
