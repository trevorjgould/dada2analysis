library(dada2)
seqtab.nochim <- readRDS("seqtab_nochim.rds")

#TAXONOMY
ref <- "/home/umii/public/dada2_taxonomy_references/sh_general_release_dynamic_s_all_25.07.2023_dev.fasta"
taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = TRUE)
saveRDS(taxa, file = "taxa.rds")

#COMBINE
both <- cbind(t(seqtab.nochim),taxa)
write.table(both, file = "ITS_combined_sequences_taxa.txt", sep = "\t", quote = FALSE, col.names=NA)