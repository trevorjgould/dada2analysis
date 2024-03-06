library(dada2)
seqtab.nochim <- readRDS("seqtab_nochim.rds")
#TAXONOMY
taxardp <- assignTaxonomy(seqtab.nochim, "/home/umii/goul0109/taxonomy/rdp_train_set_18.fa.gz", multithread=TRUE)
taxardp<- addSpecies(taxardp, "/home/umii/goul0109/taxonomy/rdp_species_assignment_18.fa.gz")
saveRDS(taxardp, file = "taxIDrdp.rds")

#TAXONOMY
taxasilva <- assignTaxonomy(seqtab.nochim, "/home/umii/goul0109/taxonomy/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxasilva <- addSpecies(taxasilva, "/home/umii/goul0109/taxonomy/silva_species_assignment_v138.1.fa")
saveRDS(taxasilva, file = "taxIDsilva.rds")

both1 <- cbind(t(seqtab.nochim),taxasilva)
write.table(both1, file = "combined_sequences_taxa_silva.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxardp)
write.table(both2, file = "combined_sequences_taxa_rdp.txt", sep = "\t", quote = FALSE, col.names=NA)