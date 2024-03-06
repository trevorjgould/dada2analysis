library(dada2)
seqtab.nochim <- readRDS("seqtab.nochim.rds")
ref <- "/home/umii/goul0109/sh_general_release_dynamic_s_all_25.07.2023_dev.fasta"
taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = TRUE)
saveRDS(taxa, file = "taxa.rds")
#
both <- cbind(t(seqtab.nochim),taxa)
write.table(both, file = "combined_sequences_taxa.txt", sep = "\t", quote = FALSE, col.names=NA)