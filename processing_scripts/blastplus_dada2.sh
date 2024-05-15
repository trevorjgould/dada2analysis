#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=60g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

cd /scratch.global/goul0109/dir/
module load ncbi_blast+
blastn -db nt -query sequences.fasta -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -num_alignments 1 > blastN_out.tab
echo -e "qseqid\sseqid\evalue\bitscore\sgi\sacc\staxids\sscinames\scomnames\stitle" > standard_blast_fields.tsv
cat standard_blast_fields.tsv blastN_out.tab > t && mv t blastN_out_tabbed.tab