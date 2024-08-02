#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=128
#SBATCH --mem=200gb
#SBATCH --tmp=200g
#SBATCH --account=kennedyp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

cd /scratch.global/goul0109/kennedy/mycophenSoil18S/

# remove adapters
module load cutadapt
module load parallel
module load R/4.4.0-openblas-rocky8

mkdir 01_adapter
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 --minimum-length 150 -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT -o 01_adapter/${i} -p 01_adapter/${i//_R1_/_R2_} ${i} ${i//_R1_/_R2_} > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" >> run_cutadapt.sh; done
chmod +x run_cutadapt.sh
parallel < run_cutadapt.sh
cd 01_adapter    
mkdir 01_logs
mv cutadapt* 01_logs
grep "passing" 01_logs/* > summary_adapter_trimming.txt

# remove primers
mkdir ../02_filtered  
#18S
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 -O 12 --minimum-length 100 --discard-untrimmed -g GTACACACCGCCCGTC -G TGATCCTTCTGCAGGTTCACCTAC -a GTAGGTGAACCTGCAGAAGGATCA -A GACGGGCGGTGTGTAC --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_/_R2_} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#18S WANDA1 AML2
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 -O 12 --minimum-length 100 --discard-untrimmed -g CAGCCGCGGTAATTCCAGCT -G GAACCCAAACACTTTGGTTTCC -a AGCTGGAATTACCGCGGCTG -A GGAAACCAAAGTGTTTGGGTTC --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_/_R2_} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
# PR2 primers
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 -O 12 --minimum-length 100 --discard-untrimmed -g TTAAARVGYTCGTAGTYG -G CCGTCAATTHCTTYAART --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_/_R2_} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
grep "passing" 02_logs/* > summary_primer_trimming.txt

Rscript ~/home/umii/goul0109/scripts/dada2_processing/dada2_18S_maar.R


for i in 2*1-S*.gz; do echo "cutadapt --cores 8 --minimum-length 100 --discard-untrimmed -g CAGCCGCGGTAATTCCAGCT -G GAACCCAAACACTTTGGTTTCC -a AGCTGGAATTACCGCGGCTG -A GGAAACCAAAGTGTTTGGGTTC --discard-untrimmed -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_/_R2_} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt18s.cmd; done
for i in 2*4-S*.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -g TCGATGAAGAACGCAGCG -G TCCTCCGCTTATTGATATGC -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done


# 
library(dada2)
path <- (".")
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=5, minLen = 100, maxEE=c(2,4), matchIDs=TRUE, maxN = 0, rm.phix=TRUE, multithread=TRUE, verbose = TRUE)
head(out)

#dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# error models
errF <- learnErrors(derep_forward, multithread=8, randomize=TRUE)
errR <- learnErrors(derep_reverse, multithread=8, randomize=TRUE)

dadaFs <- dada(derep_forward, err=errF, multithread=8, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR, multithread=8, pool="pseudo")

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, justConcatenate=TRUE)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "seqtab.rds")

# Assumes seqtab is your sequence table of merged sequences
MINLEN <- 400
MAXLEN <- 600
seqlens <- nchar(getSequences(seqtab))
seqtab.filt <- seqtab[,seqlens >= MINLEN & seqlens <= MAXLEN]

seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "seqtab_nochim.rds")
uniquesToFasta(seqtab.nochim, fout = "sequences.fasta")

  # set a little function
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
               filtered=out[,2], dada_f=sapply(dadaFs, getN),
               dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = "sequence_process_summary.txt", sep = "\t", quote=FALSE)

taxrefa <- "/home/umii/public/dada2_taxonomy_references/maarjam_dada2.txt"
taxa <- assignTaxonomy(seqtab.nochim, taxrefa, tryRC = TRUE, taxLevels = c("Class", "Order", "Family", "Genus", "Species"), multithread = TRUE, outputBootstraps = TRUE)
taxout <- taxa$tax
bootout <- taxa$boot
saveRDS(taxout, file = "taxIDmaar.rds")
saveRDS(bootout, file = "taxIDmaar_bootstrap.rds")
both1 <- cbind(t(seqtab.nochim),taxout,bootout)
write.table(both1, file = "18S_combined_sequences_taxamaar_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxout)
write.table(both2, file = "18S_combined_sequences_taxamaar.txt", sep = "\t", quote = FALSE, col.names=NA)

#library(dada2)
#seqtab.nochim <- readRDS("seqtab_nochim.rds")

#taxrefb <- "/home/umii/public/dada2_taxonomy_references/pr2_version_5.0.0_SSU_dada2.fasta.gz" 
#taxaPR2 <- assignTaxonomy(seqtab.nochim, taxrefb, multithread=TRUE, minBoot = 95, verbose = TRUE, taxLevels=c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"), outputBootstraps = TRUE)
#taxoutpr2 <- taxaPR2$tax
#bootoutpr2 <- taxaPR2$boot
#saveRDS(taxoutpr2, file = "taxIDPR2.rds")
#saveRDS(bootoutpr2, file = "taxIDPR2_bootstrap.rds")
#both1 <- cbind(t(seqtab.nochim),taxoutpr2,bootoutpr2)
#write.table(both1, file = "18S_combined_sequences_taxaPR2_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
#both2 <- cbind(t(seqtab.nochim),taxoutpr2)
#write.table(both2, file = "18S_combined_sequences_taxaPR2.txt", sep = "\t", quote = FALSE, col.names=NA)
quit("no")