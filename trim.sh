# remove adapters
mkdir 01_adapter
for i in *_R1_001.fastq.gz; do echo "~/.local/bin/cutadapt --cores 4 --minimum-length 150 -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT -o 01_adapter/${i//_R1_001.fastq.gz/-cut_R1_001.fastq.gz} -p 01_adapter/${i//_R1_001.fastq.gz/-cut_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" >> run_cutadapt.sh; done
chmod +x run_cutadapt.sh
./run_cutadapt.sh
cd 01_adapter    
mkdir 01_logs
mv cutadapt* 01_logs
#
# remove primers
mkdir ../02_filtered  
#V4
#for i in *_R1_001.fastq.gz; do echo "~/.local/bin/cutadapt --cores 4 --minimum-length 100 --discard-untrimmed -g GGACTACHVGGGTWTCTAAT -G GGACTACHVGGGTWTCTAAT --discard-untrimmed -o ../02_filtered/${i//-cut_R1_001.fastq.gz/-trimmed_R1_001.fastq.gz} -p ../02_filtered/${i//-cut_R1_001.fastq.gz/-trimmed_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#V3V4
for i in *_R1_001.fastq.gz; do echo "~/.local/bin/cutadapt --cores 4 --minimum-length 100 --discard-untrimmed -g AGAGTTTGATCMTGGCTCAG -G ATTACCGCGGCTGCTGG --discard-untrimmed -o ../02_filtered/${i//-cut_R1_001.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//-cut_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
chmod +x run_cutadapt2.cmd
./run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
