#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=8
#SBATCH --mem=60
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

cd /location/of/fastq.gz/files/
module load R/4.1.0
./trim_18S.sh
cd 02_filtered/
Rscript ../run_18S_dada2.R