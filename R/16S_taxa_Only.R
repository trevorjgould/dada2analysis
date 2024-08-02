library(dada2)
seqtab.nochim <- readRDS("seqtab_nochim.rds")

#TAXONOMY
taxasilva <- assignTaxonomy(seqtab.nochim, "/home/umii/public/dada2_taxonomy_references/silva_nr99_v138.1_train_set.fa", multithread=TRUE, outputBootstraps = TRUE)
taxasilva$tax <- addSpecies(taxasilva$tax, "/home/umii/public/dada2_taxonomy_references/silva_species_assignment_v138.1.fa")
saveRDS(taxasilva$tax, file = "taxIDsilva.rds")
saveRDS(taxasilva$boot, file = "taxIDsilva_bootstrap.rds")
both1 <- cbind(t(seqtab.nochim),taxasilva$tax, taxasilva$boot)
both2 <- cbind(t(seqtab.nochim),taxasilva$tax)
write.table(both1, file = "16S_combined_sequences_taxa_silva_boot.txt", sep = "\t", quote = FALSE, col.names=NA)
write.table(both2, file = "16S_combined_sequences_taxa_silva.txt", sep = "\t", quote = FALSE, col.names=NA)
quit("no")