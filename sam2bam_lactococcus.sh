#!/bin/bash -l
#SBATCH -D /home/smhigdon/Projects/lactococcus
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/sam2bam_lactococcus-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/sam2bam_lactococcus-stderr-%j.txt
#SBATCH -J sam2bam_lacto
#SBATCH --array 3-72 # number of files to be processed
#SBATCH -p high
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=4000
#SBATCH -t 96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu

# Name: sam2bam_lactococcus.sh
# Created by: Shawn Higdon
# creation date: June 23, 2020
# Description: An array script to convert sam files to bam files.

# Load modules
module load samtools

# Define Variable

OUTPUT_FOLDER=analysis # folder for all output
mkdir -p $OUTPUT_FOLDER

SEEDFILE=input_files/sam2bam_inputs.txt # the list of .sam files that will be processed (absolute path), 1 per line.

SEED=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

bname=`basename $SEED`
bname=`echo $bname | cut -d. -f1`

SAMPLE_FOLDER=$OUTPUT_FOLDER/$bname
mkdir -p $SAMPLE_FOLDER

MAPPING_FOLDER=$SAMPLE_FOLDER/bwa_mapping
mkdir -p $MAPPING_FOLDER

ASM=$SAMPLE_FOLER/megahit/${bname}.contigs.fa

# Step One: index the assembly with samtools

samtools faidx $ASM

# Step Two: convert the .sam file to .bam file

samtools import $ASM $SEED $MAPPING_FOLDER/${bname}.bam

samtools sort $MAPPING_FOLDER/${bname}.bam -o $MAPPING_FOLDER/${bname}.bam.sorted

samtools index $MAPPING_FOLDER/${bname}.bam.sorted

hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks
stream
