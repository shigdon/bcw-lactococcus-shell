#!/bin/bash -l
#SBATCH -D /home/smhigdon/Projects/lactococcus
#SBATCH -o /home/smhigdon/Projects/lactococcus/slurm-log/sm_compute-stdout-%j.txt
#SBATCH -e /home/smhigdon/Projects/lactococcus/slurm-log/sm_compute-stderr-%j.txt
#SBATCH -J sm_compute
#SBATCH --array 1-81 # number of files to be processed
#SBATCH -p high
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=8000
#SBATCH -t 96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smhigdon@ucdavis.edu

# Name: sm_compute.sh
# Created by: Shawn Higdon
# creation date: Sep 6, 2019
# Description: A script to compute Minhash signatures from MEGAhit assemblies of trimmed PE150 Illumina reads that correspond to isolate genomes using Sourmash - as an array job.

# Load modules

source activate sourmash

# Define Variable

OUTPUT_FOLDER=analysis # folder for all output
mkdir -p $OUTPUT_FOLDER

SEEDFILE=input_files/sourmash_compute_inputs.txt # the list of assembly files that will be processed (absolute path), 1 per line.

SEED=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

bname=`basename $SEED`
bname=`echo $bname | cut -d. -f1`

SAMPLE_FOLDER=$OUTPUT_FOLDER/${bname}
mkdir -p $SAMPLE_FOLDER

SOURMASH_FOLDER=$SAMPLE_FOLDER/sourmash
mkdir -p $SOURMASH_FOLDER

#  compute Minhash Sketches with sourmash v. 2.0.11

sourmash compute -k 21 --scaled 2000 -o $SOURMASH_FOLDER/${bname}.k21.sig -f $SEED

sourmash compute -k 31 --scaled 2000 -o $SOURMASH_FOLDER/${bname}.k31.sig -f $SEED

sourmash compute -k 51 --scaled 2000 -o $SOURMASH_FOLDER/${bname}.k51.sig -f $SEED

source deactivate

hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks
stream
