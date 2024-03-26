#!/bin/bash -l        
#SBATCH --time=12:00:00
#SBATCH --ntasks=8
#SBATCH --mem=60gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=goul0109@umn.edu

cd /location/of/fastq.gz/files/
module load R/4.1.0
./trim_16S.sh
cd 02_filtered/
Rscript ../run_16S_dada2.R