# illumina index seq primer
GTAGGTGAACCTGCAGAAGGATCA TG CTGACTGACT

# 18S euk R1 primer
GTACACACCGCCCGTC
# 18S euk R2 primer
TGATCCTTCTGCAGGTTCACCTAC
/home/umii/goul0109/scripts/primer_check2.sh {names of R1 file here}
/home/umii/goul0109/scripts/fix_names.sh *
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

module load fastqc 
fastqc a few files to understand length and quality

module load R
R

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

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=10)

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


  # set a little function
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
               filtered=out[,2], dada_f=sapply(dadaFs, getN),
               dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = "sequence_process_summary.txt", sep = "\t", quote=FALSE)

taxrefa <- "/home/umii/goul0109/maarjam_dada2.txt"
taxa <- assignTaxonomy(seqtab.nochim, taxrefa, tryRC = TRUE, taxLevels = c("Class", "Order", "Family", "Genus", "Species"), multithread = TRUE, outputBootstraps = TRUE)
saveRDS(taxa, file = "18SmaarjamtaxID.rds")

#taxrefb <- "/home/kennedyp/shared/taxonomy/pr2_version_5.0.0_SSU_dada2.fasta.gz" 
#taxaPR2 <- assignTaxonomy(seqtab.nochim, taxrefb, multithread=TRUE, minBoot = 95, verbose = TRUE, taxLevels=c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"), outputBootstraps = TRUE)
#saveRDS(taxaPR2, file = "18StaxID.rds")

both <- cbind(t(seqtab.nochim),taxa)
write.table(both, file = "combined_sequences_taxa.txt", sep = "\t", quote = FALSE, col.names=NA)