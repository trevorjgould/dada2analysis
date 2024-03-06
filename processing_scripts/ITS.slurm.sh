#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

cd /scratch.global/goul0109/kennedy/B4W_Necromass/ITS/
module load R/4.1.0
./trim_ITS.sh
cd /scratch.global/goul0109/kennedy/B4W_Necromass/ITS/02_filtered/
Rscript run_ITS_dada2.R
