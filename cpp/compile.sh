#!/bin/bash
#SBATCH -c 1  # Number of Cores per Task
#SBATCH -p cpu # Partition
#SBATCH -t 00:10:00 # Job time limit
#SBATCH --mem=8G  # Requested Memory
#SBATCH -o slurm-%j.sout  # %j = job ID
#SBATCH --constraint=avx512 

module load uri/main
module load GSL/2.7-GCC-11.2.0

make
