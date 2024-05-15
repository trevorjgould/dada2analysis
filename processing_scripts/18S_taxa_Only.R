library(dada2)
seqtab.nochim <- readRDS("seqtab_nochim.rds")

taxrefb <- "/home/umii/public/dada2_taxonomy_references/pr2_version_5.0.0_SSU_dada2.fasta.gz" 
taxaPR2 <- assignTaxonomy(seqtab.nochim, taxrefb, multithread=TRUE, minBoot = 95, verbose = TRUE, taxLevels=c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"), outputBootstraps = TRUE)

taxout <- taxaPR2$tax
bootout <- taxaPR2$boot
saveRDS(taxout, file = "taxIDPR2.rds")
saveRDS(bootout, file = "taxIDPR2_bootstrap.rds")

#saveRDS(taxaPR2$tax, file = "18StaxID.rds")
#saveRDS(taxaPR2$boot, file = "18StaxID_bootstrap.rds")
both1 <- cbind(t(seqtab.nochim),taxaPR2$tax,taxaPR2$boot)
write.table(both1, file = "18S_combined_sequences_taxaPR2_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxaPR2$tax)
write.table(both2, file = "18S_combined_sequences_taxaPR2.txt", sep = "\t", quote = FALSE, col.names=NA)
quit("no")

