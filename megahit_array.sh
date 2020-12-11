#!/bin/bash -l
#SBATCH -D /home/smhigdon/Projects/lactococcus
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/megahit_array-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/megahit_array-stderr-%j.txt
#SBATCH -J megahit
#SBATCH --array 1-81 # number of files to be processed
#SBATCH -p high
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8000
#SBATCH -t 96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu

# Name: megahit_array.sh
# Created by: Shawn Higdon
# creation date: Sep 05, 2019
# Description: A script to assemble trimmed PE150 Illumina reads for isolate genomes with megahit and assess assembly quality using QUAST - as an array job.

# Load modules

module load bio

# Define Variable

OUTPUT_FOLDER=analysis # folder for all output
mkdir -p $OUTPUT_FOLDER

SEEDFILE=input_files/megahit_R1.txt # megahit_bcw_unknowns_input_1.txt is the list of PE_1 files that will be processed (absolute path), 1 per line.

SEED=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

SEEDFILE2=input_files/megahit_R2.txt # megahit_bcw_unknowns_input_2.txt is the list of PE_2 files that will be processed (absolute path), 1 per line.

SEED2=$(cat $SEEDFILE2 | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

bname=`basename $SEED`
bname=`echo $bname | cut -d. -f1`

SAMPLE_FOLDER=$OUTPUT_FOLDER/${bname}
mkdir -p $SAMPLE_FOLDER

#  assemble the trimmed read files with megahit v 1.1.1

MEGAHIT_FOLDER=$SAMPLE_FOLDER/megahit
megahit -1 $SEED -2 $SEED2 --out-dir $MEGAHIT_FOLDER --out-prefix ${bname} -t 4 

# assess genome assembly quality with quast

QUAST_FOLDER=$SAMPLE_FOLDER/quast
quast $MEGAHIT_FOLDER/${bname}.contigs.fa --output-dir $QUAST_FOLDER --labels ${bname}

hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks
stream
