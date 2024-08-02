#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=128
#SBATCH --mem=499gb
#SBATCH --tmp=200g
#SBATCH --account=kennedyp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

cd /scratch.global/goul0109/kennedy/amazITS/

# remove adapters
module load cutadapt
module load parallel
module load R/4.4.0-openblas-rocky8

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
#ITS2
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 -O 12 --pair-filter=any --minimum-length 100 --discard-untrimmed -g TCGATGAAGAACGCAGCG -G TCCTCCGCTTATTGATATGC -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#ITS1
#for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 -O 12 --pair-filter=any --minimum-length 100 --discard-untrimmed -g CTTGGTCATTTAGAGGAAGTAA -G GCTGCGTTCTTCATCGATGC -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
#ITScustom
for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 -O 12 --pair-filter=any --minimum-length 100 --discard-untrimmed  -g AACTTTYRRCAAYGGATCWCT -G AGCCTCCGCTTATTGATATGCTTAART -o ../02_filtered/${i} -p ../02_filtered/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done

chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
grep "passing" 02_logs/* > summary_primer_trimming.txt

Rscript ~/home/umii/goul0109/scripts/dada2_processing/DADA2_ITS_128.R
####################################################

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

# ITSs
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# poor quality
out <- filterAndTrim(fnFs[1:2], filtFs[1:2], fnRs[1:2], filtRs[1:2], maxN = 0, maxEE = c(4, 6), truncQ = 2, minLen = 100, truncLen=c(250,175), rm.phix = TRUE, compress = TRUE, multithread = TRUE)
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

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=20)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "seqtab.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=8, verbose=TRUE)
dim(seqtab.nochim)
lname <- nchar(colnames(seqtab.nochim))
seqtab.nochim <- seqtab.nochim[,(lname > 280)]
saveRDS(seqtab.nochim, "seqtab_nochim.rds")

uniquesToFasta(seqtab.nochim, fout = "sequences.fasta")

  # set a little function
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
               filtered=out[,2], dada_f=sapply(dadaFs, getN),
               dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = "sequence_process_summary.txt", sep = "\t", quote=FALSE)

ref <- "/home/umii/public/dada2_taxonomy_references/sh_general_release_dynamic_all_04.04.2024.fasta"
taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = TRUE, outputBootstraps = TRUE)

taxout <- taxa$tax
bootout <- taxa$boot
saveRDS(taxout, file = "taxID.rds")
saveRDS(bootout, file = "taxID_bootstrap.rds")

#saveRDS(taxa$tax, file = "taxa.rds")
#saveRDS(taxa$boot, file = "taxa_bootstrap.rds")
#
both1 <- cbind(t(seqtab.nochim),taxa$tax, taxa$boot)
write.table(both1, file = "ITS_combined_sequences_taxa_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxa$tax)
write.table(both2, file = "ITS_combined_sequences_taxa.txt", sep = "\t", quote = FALSE, col.names=NA)

####################################################

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=60gb
#SBATCH --tmp=200g
#SBATCH --account=kennedyp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

cd /scratch.global/goul0109/kennedy/AmazonSoilITS2/02_filtered/
module load R
Rscript ITS_taxa.R

library(dada2)
seqtab.nochim <- readRDS("seqtab_nochim.rds")

#TAXONOMY
ref <- "/home/umii/public/dada2_taxonomy_references/sh_general_release_dynamic_all_04.04.2024.fasta"
taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = 8, outputBootstraps = TRUE)

taxout <- taxa$tax
bootout <- taxa$boot
saveRDS(taxout, file = "taxID.rds")
saveRDS(bootout, file = "taxID_bootstrap.rds")

#saveRDS(taxa$tax, file = "taxa.rds")
#saveRDS(taxa$boot, file = "taxa_bootstrap.rds")
#
both1 <- cbind(t(seqtab.nochim),taxa$tax, taxa$boot)
write.table(both1, file = "ITS_combined_sequences_taxa_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxa$tax)
write.table(both2, file = "ITS_combined_sequences_taxa.txt", sep = "\t", quote = FALSE, col.names=NA)

####################################################

mkdir ../DADA2_OUTPUT
mv *.rds ../DADA2_OUTPUT/
mv *.txt ../DADA2_OUTPUT/
cd ..
mv DADA2_OUTPUT/ XXXX_DADA2_OUTPUT/