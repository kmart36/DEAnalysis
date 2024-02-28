#!/bin/bash 
#SBATCH -p medium # partition (queue) 
#SBATCH -N 1 # (leave at 1 unless using multi-node specific code) 
#SBATCH -n 1 # number of cores 
#SBATCH --mem-per-cpu=16384 # memory per core 
#SBATCH --job-name="nextflow" # job name 
#SBATCH -o slurm.%N.%j.stdout.txt # STDOUT 
#SBATCH -e slurm.%N.%j.stderr.txt # STDERR 
#SBATCH --mail-user=kam071@bucknell.edu # address to email 
#SBATCH --mail-type=ALL # mail events (NONE, BEGIN, END, FAIL, ALL) 

module load qc
module load kraken2
module load assembly_assessment
module load ncbi-blast
conda activate cd-hit
/home/kam071/nextflow run DEanalysis.nf -with-report -with-timeline -with-dag
