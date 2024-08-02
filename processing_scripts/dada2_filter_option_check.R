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

# 16s
out1reg <- filterAndTrim(fnFs[1:2], filtFs[1:2], fnRs[1:2], filtRs[1:2], maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
out2ee24 <- filterAndTrim(fnFs[1:2], filtFs[1:2], fnRs[1:2], filtRs[1:2], maxN = 0, maxEE = c(2, 4), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
out2eeinf <- filterAndTrim(fnFs[1:2], filtFs[1:2], fnRs[1:2], filtRs[1:2], maxN = 0, truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
out1nomin <- filterAndTrim(fnFs[1:2], filtFs[1:2], fnRs[1:2], filtRs[1:2], maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
out1nophix <- filterAndTrim(fnFs[1:2], filtFs[1:2], fnRs[1:2], filtRs[1:2], maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, compress = TRUE, multithread = TRUE)
out1nomaxN <- filterAndTrim(fnFs[1:2], filtFs[1:2], fnRs[1:2], filtRs[1:2], maxEE = c(2, 2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out1reg)
head(out2ee24)
head(out2eeinf)
head(out1nomin)
head(out1nophix)
head(out1nomaxN)