library(dada2)
seqtab.nochim <- readRDS("seqtab_nochim.rds")

#TAXONOMY
ref <- "/home/umii/public/dada2_taxonomy_references/sh_general_release_dynamic_s_all_25.07.2023_dev.fasta"
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
