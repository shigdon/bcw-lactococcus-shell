#!/bin/bash -l

#SBATCH -p high
#SBATCH -D /home/smhigdon/Projects/lactococcus
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/trim-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/trim-stderr-%j.txt
#SBATCH -J bcw_trim
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --array 76 
#SBATCH -t 96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu

# Name: trimmomatic_array.sh
# Created by: Shawn Higdon
# creation date: Sep 05, 2019
# A script to quality trim and debarcode fastq reads.

# Load modules

module load bio

# Define Variable

OUTPUT_FOLDER=analysis # folder for all output
mkdir -p $OUTPUT_FOLDER

SEEDFILE=input_files/lactococcus_isolates_R1_input.txt # the list of PE_1 files that will be processed (absolute path), 1 per line.

SEED=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

SEEDFILE2=input_files/lactococcus_isolates_R2_input.txt # the list of PE_2 files that will be processed (absolute path), 1 per line.

SEED2=$(cat $SEEDFILE2 | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

bname=`basename $SEED`
bname=`echo $bname | cut -d_ -f1`

SAMPLE_FOLDER=$OUTPUT_FOLDER/${bname}
mkdir -p $SAMPLE_FOLDER

# Trim the read files with Trimmomatic v0.36

trimmomatic PE \
	-threads 12 \
	$SEED \
	$SEED2 \
	$SAMPLE_FOLDER/${bname}.trim_R1.fq.gz \
	$SAMPLE_FOLDER/${bname}.orphan_R1.fq.gz \
	$SAMPLE_FOLDER/${bname}.trim_R2.fq.gz \
	$SAMPLE_FOLDER/${bname}.orphan_R2.fq.gz \
	ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15:8:TRUE \
	LEADING:2 \
	TRAILING:2 \
	SLIDINGWINDOW:4:15 \
	MINLEN:50

hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks
stream
