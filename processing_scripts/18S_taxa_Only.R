library(dada2)
seqtab.nochim <- readRDS("seqtab_nochim.rds")

taxrefa <- "/home/umii/public/dada2_taxonomy_references/maarjam_dada2.txt"
taxa <- assignTaxonomy(seqtab.nochim, taxrefa, tryRC = TRUE, taxLevels = c("Class", "Order", "Family", "Genus", "Species"), multithread = TRUE, outputBootstraps = TRUE)
taxout <- taxa$tax
bootout <- taxa$boot
saveRDS(taxout, file = "taxIDmaar.rds")
saveRDS(bootout, file = "taxIDmaar_bootstrap.rds")
both1 <- cbind(t(seqtab.nochim),taxout,bootout)
write.table(both1, file = "18S_combined_sequences_taxamaar_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxout)
write.table(both2, file = "18S_combined_sequences_taxamaar.txt", sep = "\t", quote = FALSE, col.names=NA)

taxrefb <- "/home/umii/public/dada2_taxonomy_references/pr2_version_5.0.0_SSU_dada2.fasta.gz" 
taxaPR2 <- assignTaxonomy(seqtab.nochim, taxrefb, multithread=TRUE, minBoot = 95, verbose = TRUE, taxLevels=c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"), outputBootstraps = TRUE)
taxoutpr2 <- taxaPR2$tax
bootoutpr2 <- taxaPR2$boot
saveRDS(taxoutpr2, file = "taxIDPR2.rds")
saveRDS(bootoutpr2, file = "taxIDPR2_bootstrap.rds")
both1 <- cbind(t(seqtab.nochim),taxoutpr2,bootoutpr2)
write.table(both1, file = "18S_combined_sequences_taxaPR2_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxoutpr2)
write.table(both2, file = "18S_combined_sequences_taxaPR2.txt", sep = "\t", quote = FALSE, col.names=NA)
quit("no")