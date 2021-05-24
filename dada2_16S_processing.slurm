#!/bin/bash -l        
#SBATCH --time=6:00:00
#SBATCH --ntasks=12
#SBATCH --mem=30gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=goul0109@umn.edu

cd /panfs/roc/scratch/goul0109/USER/16S/02_filtered/
module load R/3.6.0
Rscript /panfs/roc/scratch/goul0109/USER/dada2_pipeline/trim.sh

Rscript /panfs/roc/scratch/goul0109/USER/dada2_pipeline/dada2_Only.R
