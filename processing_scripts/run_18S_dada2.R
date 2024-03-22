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

taxrefa <- "/home/umii/public/dada2_taxonomy_references/maarjam_dada2.txt"
taxa <- assignTaxonomy(seqtab.nochim, taxrefa, tryRC = TRUE, taxLevels = c("Class", "Order", "Family", "Genus", "Species"), multithread = TRUE, outputBootstraps = TRUE)
saveRDS(taxa, file = "18SmaarjamtaxID.rds")

taxrefb <- "/home/umii/public/dada2_taxonomy_references/pr2_version_5.0.0_SSU_dada2.fasta.gz" 
taxaPR2 <- assignTaxonomy(seqtab.nochim, taxrefb, multithread=TRUE, minBoot = 95, verbose = TRUE, taxLevels=c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"), outputBootstraps = TRUE)
saveRDS(taxaPR2, file = "18StaxID.rds")

both <- cbind(t(seqtab.nochim),taxa)
write.table(both, file = "18S_combined_sequences_taxamaarjam.txt", sep = "\t", quote = FALSE, col.names=NA)

both2 <- cbind(t(seqtab.nochim),taxaPR2)
write.table(both2, file = "18S_combined_sequences_taxaPR2.txt", sep = "\t", quote = FALSE, col.names=NA)
