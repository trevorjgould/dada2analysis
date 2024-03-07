##########################
# Data processing steps 
# remove adapters
module load cutadapt
mkdir 01_adapter
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 4 --minimum-length 150 -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT -o 01_adapter/${i} -p 01_adapter/${i//_R1_/_R2_} ${i} ${i//_R1_/_R2_} > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" >> run_cutadapt.sh; done
chmod +x run_cutadapt.sh
./run_cutadapt.sh
cd 01_adapter    
mkdir 01_logs
mv cutadapt* 01_logs
grep "passing" 01_logs/* > summary_adapter_trimming.txt


# remove primers
mkdir ../02_filtered  
#18S
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 4 --minimum-length 100 --discard-untrimmed -g GTACACACCGCCCGTC -G TGATCCTTCTGCAGGTTCACCTAC -a GTAGGTGAACCTGCAGAAGGATCA -A GACGGGCGGTGTGTAC --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_/_R2_} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#18S WANDA1 AML2
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 4 --minimum-length 100 --discard-untrimmed -g CAGCCGCGGTAATTCCAGCT -G GAACCCAAACACTTTGGTTTCC -a AGCTGGAATTACCGCGGCTG -A GGAAACCAAAGTGTTTGGGTTC --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_/_R2_} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done


chmod +x run_cutadapt2.cmd
./run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
grep "passing" 02_logs/* > summary_primer_trimming.txt
