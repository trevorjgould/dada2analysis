#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=499g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

cd /home/umii/goul0109/scratch.global/goul0109/kennedy/Final/Corn/02_filtered/
module load ncbi_blast+
blastn -db nt -query sequences.fasta -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle'  -max_target_seqs 1  -max_hsps 1 -num_threads 128 > blastN_out.tab
cat /home/umii/goul0109/scripts/standard_blast_fields.tsv blastN_out.tab > t && mv t blastN_out_tabbed.tab
