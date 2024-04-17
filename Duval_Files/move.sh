#!/bin/bash 
#SBATCH -p short # partition (queue) 
#SBATCH -N 1 # (leave at 1 unless using multi-node specific code) 
#SBATCH -n 1 # number of cores 
#SBATCH --mem-per-cpu=8192 # memory per core 
#SBATCH --job-name="nextflow" # job name 
#SBATCH -o slurm.%N.%j.stdout.txt # STDOUT 
#SBATCH -e slurm.%N.%j.stderr.txt # STDERR 
#SBATCH --mail-user=kam071@bucknell.edu # address to email 
#SBATCH --mail-type=NONE # mail events (NONE, BEGIN, END, FAIL, ALL) 

mv fastqc_files/*/*.zip fastqc_summaries