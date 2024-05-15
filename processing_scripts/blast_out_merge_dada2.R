# blast output header is just V1...V10

# V1 has sequence number and more than one occasionally skipping a number
blastN_out_tabbed <- read.delim("blastN_out_tabbed.tab", header=FALSE)
# now just one per line still skipping missing
blastNshort <- blastN_out_tabbed[!duplicated(blastN_out_tabbed$V1),]

# add row number column to dada2 output
seqtab <- readRDS("seqtab_nochim.rds")
seqtab <- t(seqtab)
seqtab <- as.data.frame(seqtab)
seqtab$V1 <- 1:nrow(seqtab)

# merge by V1
both <- merge(seqtab,blastNshort, by = "V1", all.x=TRUE)
rownames(both) <- rownames(seqtab)
View(both)

taxa <- readRDS("taxIDsilva.rds")
both2 <- merge(both,taxa, by = "row.names")
both2 <- both2[order(both2$V1),]
both2 <- both2[,-c(2)]
write.table(both2, file = "dada2_blast_silva.txt", sep = "\t", quote = FALSE)
