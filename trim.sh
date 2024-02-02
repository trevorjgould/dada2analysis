#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=64
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

cd /scratch.global/goul0109/

# remove adapters
module load cutadapt
module load parallel
mkdir 01_adapter
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 150 -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT -o 01_adapter/${i//_R1_001.fastq.gz/_R1_001.fastq.gz} -p 01_adapter/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" >> run_cutadapt.sh; done
chmod +x run_cutadapt.sh
parallel < run_cutadapt.sh
cd 01_adapter    
mkdir 01_logs
mv cutadapt* 01_logs
grep "passing" 01_logs/* > summary_adapter_trimming.txt
#
# remove primers
mkdir ../02_filtered  
#V4
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -g GTGCCAGCMGCCGCGGTAA -G SGACTACHVGGGTWTCTAAT --discard-untrimmed -o ../02_filtered/${i//_R1_001.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#V3V4
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -g CCTACGGGNGGCWGCAG -G NGACTACHVGGGTWTCTAAT --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#ITS
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 4 --pair-filter=any --minimum-length 100 --discard-untrimmed -g CTTGGTCATTTAGAGGAAGTAA -G GCTGCGTTCTTCATCGATGC --discard-untrimmed -o ../02_filtered/${i//-cut_R1_001.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//-cut_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#V1V3
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -g AGAGTTTGATCMTGGCTCAG -G ATTACCGCGGCTGCTGG --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done

chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
grep "passing" 02_logs/* > summary_primer_trimming.txt
#

for i in Oral-Rinse-16S-rRNA-V1V3-S103-L001mapped_and_unmapped_bothReadsUnmapped_host_removed_R1.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -g AGAGTTTGATCMTGGCTCAG -G ATTACCGCGGCTGCTGG --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//R1.fastq.gz/_R2.fastq.gz} ${i} ${i//R1_/R2_} > ../02_filtered/cutadapt.${i//R1.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
for i in Oral-Rinse-16S-rRNA-V3V4-S106-L001mapped_and_unmapped_bothReadsUnmapped_host_removed_R1.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -g CCTACGGGNGGCWGCAG -G NGACTACHVGGGTWTCTAAT --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//R1.fastq.gz/_R2.fastq.gz} ${i} ${i//R1_/R2_} > ../02_filtered/cutadapt.${i//R1.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done

library(dada2)
path <- (".")
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# 16s
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200), maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=8)
out


# removing nextera transposase adapter contamination

>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

mkdir 03_NTadapter
for i in *_R1.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 150 -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -G CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o 03_NTadapter/${i} -p 03_NTadapter/${i//_R1/_R2} ${i} ${i//_R1/_R2} > 03_NTadapter/cutadapt.${i//_R1.fastq.gz/.log.txt}" >> run_cutadaptNT.sh; done
chmod +x run_cutadaptNT.sh
parallel < run_cutadaptNT.sh
cd 03_NTadapter    
mkdir 03_logs
mv cutadapt* 03_logs
grep "passing" 03_logs/* > summary_NTadapter_trimming.txt

