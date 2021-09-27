#' Create_Tables
#'
#' This function processes dada2 output and returns tables for plots and stats
#' @param inputtable table
#' @param metadata table
#' @param taxa table
#'
#' @export
#'
#' @importFrom utils write.table

Create_Tables <- function(inputtable,metadata,taxa){
#row.names(inputtable) <- gsub("-trimmed","", row.names(inputtable))
#find common names
common <- intersect(rownames(metadata),rownames(inputtable))
# get just the overlapping samples
newmap <- metadata[common,, drop = FALSE]
newtable <- inputtable[common,, drop = FALSE]
#save to file
saveRDS(newtable, file = "Sequence_table_common.rds")
write.table(newmap, file = "Metadata_common.txt")
#Taxa table
newtable1 <- t(newtable)
#find common names
common <- intersect(rownames(taxa),rownames(newtable1))
# get just the overlapping samples
taxa2 <- taxa[common,, drop = FALSE]
newtable2 <- newtable1[common,, drop = FALSE]
both <- cbind(newtable2,taxa2)
both <- as.data.frame(both)

# replacing NA with Unknown higher level taxa
both$Phylum <- ifelse(is.na(both$Phylum), paste0("Unknown ",both$Kingdom), both$Phylum)
both$Class <- ifelse(is.na(both$Class), paste0("Unknown ",both$Phylum), both$Class)
both$Order <- ifelse(is.na(both$Order), paste0("Unknown ",both$Class), both$Order)
both$Family <- ifelse(is.na(both$Family), paste0("Unknown ",both$Order), both$Family)
both$Genus <- ifelse(is.na(both$Genus), paste0("Unknown ",both$Family), both$Genus)
both$Species <- ifelse(is.na(both$Species), paste0("Unknown ",both$Genus), both$Species)
both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
row.names(both) <- row.names(newtable1)
#save to file
write.table(both, file = "combined_sequences_taxa.txt", sep = "\t", quote = FALSE)
return(list(newtable = newtable, newmap = newmap, combined_taxa = both))
}
