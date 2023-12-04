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
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 150 -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT -o 01_adapter/${i//_R1_001.fastq.gz/-cut_R1_001.fastq.gz} -p 01_adapter/${i//_R1_001.fastq.gz/-cut_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" >> run_cutadapt.sh; done
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
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 4 --pair-filter=any --minimum-length 100 --discard-untrimmed -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT --discard-untrimmed -o ../02_filtered/${i//-cut_R1_001.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//-cut_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#V3V4
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 4 --pair-filter=any --minimum-length 100 --discard-untrimmed -g CCTACGGGNGGCWGCAG -G GGACTACHVGGGTWTCTAAT --discard-untrimmed -o ../02_filtered/${i//-cut_R1_001.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//-cut_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#ITS
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 4 --pair-filter=any --minimum-length 100 --discard-untrimmed -g CTTGGTCATTTAGAGGAAGTAA -G GCTGCGTTCTTCATCGATGC --discard-untrimmed -o ../02_filtered/${i//-cut_R1_001.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//-cut_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done

chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
grep "passing" 02_logs/* > summary_primer_trimming.txt
#

