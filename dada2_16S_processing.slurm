#!/bin/bash -l        
#SBATCH --time=12:00:00
#SBATCH --ntasks=24
#SBATCH --mem=60gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=goul0109@umn.edu

cd /panfs/roc/scratch/goul0109/pragman2_Project_007/16S/02_filtered/
module load R/3.6.0

Rscript /panfs/roc/scratch/goul0109/pragman2_Project_007/dada2_pipeline/dada2_Only.R
 
