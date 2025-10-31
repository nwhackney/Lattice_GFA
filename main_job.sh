#!/bin/bash

#SBATCH -c 1
#SBATCH -p cpu-long
#SBATCH -t 36:00:00
#SBATCH --mem=8G
#SBATCH -o slurm-%j.sout
#SBATCH --job-name L250
#SBATCH --constraint=avx512 

module load uri/main
module load GSL/2.7-GCC-11.2.0

ulimit -s hard
./main