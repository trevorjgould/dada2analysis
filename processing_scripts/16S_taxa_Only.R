library(dada2)
seqtab.nochim <- readRDS("seqtab_nochim.rds")

#TAXONOMY
taxasilva <- assignTaxonomy(seqtab.nochim, "/home/umii/public/dada2_taxonomy_references/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxasilva <- addSpecies(taxasilva, "/home/umii/public/dada2_taxonomy_references/silva_species_assignment_v138.1.fa")
saveRDS(taxasilva, file = "taxIDsilva.rds")

#COMBINE
both1 <- cbind(t(seqtab.nochim),taxasilva)
write.table(both1, file = "16S_combined_sequences_taxa_silva.txt", sep = "\t", quote = FALSE, col.names=NA)