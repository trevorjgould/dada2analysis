#!/bin/bash -l        
#SBATCH --time=12:00:00
#SBATCH --ntasks=64
#SBATCH --mem=60gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=goul0109@umn.edu

cd /location/of/fastq.gz/files/
compat-exec ./trim_16S.sh
cd 02_filtered/

module load R/4.3.0-openblas-rocky8
Rscript ../run_16S_dada2.R