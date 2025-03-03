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

# check if sampleID starts with a number and if so add an S to the start.
#if(is.numeric(substring(common[[1]], 1, 1))){
#rownames(newmap) <- paste0("S",common)
#rownames(newtable) <- paste0("S",common)
#}
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

  c1 <- colnames(both)
  r1 <- rownames(both)

# replacing NA with Unknown higher level taxa
both$Phylum <- ifelse(is.na(both$Phylum), paste0("Unknown ",both$Kingdom), both$Phylum)
both$Class <- ifelse(is.na(both$Class), paste0("Unknown ",both$Phylum), both$Class)
both$Order <- ifelse(is.na(both$Order), paste0("Unknown ",both$Class), both$Order)
both$Family <- ifelse(is.na(both$Family), paste0("Unknown ",both$Order), both$Family)
both$Genus <- ifelse(is.na(both$Genus), paste0("Unknown ",both$Family), both$Genus)
if("Species" %in% colnames(both)){
  both$Species <- ifelse(is.na(both$Species), paste0("Unknown ",both$Genus), both$Species);
}
both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))

  rownames(both) <- r1
  colnames(both) <- c1
  
#save to file
write.table(both, file = "combined_sequences_taxa.txt", sep = "\t", quote = FALSE)
# check that sizes of input and output are the same
# does inputtable = newtable
print("Checking if dimensions are the same")
print("inputtable Rows")
if(dim(inputtable)[1] == dim(newtable2)[1]){print("good")}
print("inputtable Columns")
if(dim(inputtable)[2] == dim(newtable2)[2]){print("good")}
# does metadata = newmap
print("metadatarows")
if(dim(metadata)[1] == dim(newmap)[1]){print("good")}
print("metadata columns")
if(dim(metadata)[2] == dim(newmap)[2]){print("good")}
# does taxa + inputtable = combined_taxa

return(list(newtable = newtable, newmap = newmap, combined_taxa = both))
}


Create_Tables_18S <- function(inputtable,metadata,taxa){
  #row.names(inputtable) <- gsub("-trimmed","", row.names(inputtable))
  #find common names
  common <- intersect(rownames(metadata),rownames(inputtable))
  # get just the overlapping samples
  newmap <- metadata[common,, drop = FALSE]
  newtable <- inputtable[common,, drop = FALSE]
  
  # check if sampleID starts with a number and if so add an S to the start.
  #if(is.numeric(substring(common[[1]], 1, 1))){
  rownames(newmap) <- paste0("S",common)
  rownames(newtable) <- paste0("S",common)
  #}
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
  
  c1 <- colnames(both)
  r1 <- rownames(both)

 # replacing NA with Unknown higher level taxa
  both$Domain <- ifelse(is.na(both$Domain), paste0("Unknown ",both$Domain), both$Domain)
  both$Supergroup <- ifelse(is.na(both$Supergroup), paste0("Unknown ",both$Domain), both$Supergroup)
  both$Division <- ifelse(is.na(both$Division), paste0("Unknown ",both$Supergroup), both$Division)
  both$Subdivision <- ifelse(is.na(both$Subdivision), paste0("Unknown ",both$Division), both$Subdivision)
  both$Class <- ifelse(is.na(both$Class), paste0("Unknown ",both$Subdivision), both$Class)
  both$Order <- ifelse(is.na(both$Order), paste0("Unknown ",both$Class), both$Order)
  both$Family <- ifelse(is.na(both$Family), paste0("Unknown ",both$Order), both$Family)
  both$Genus <- ifelse(is.na(both$Genus), paste0("Unknown ",both$Family), both$Genus)
  both$Species <- ifelse(is.na(both$Species), paste0("Unknown ",both$Genus), both$Species)
  both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
  both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
  both <- data.frame(lapply(both, function(x) {gsub("Unknown Unknown", "Unknown", x)}))

  rownames(both) <- r1
  colnames(both) <- c1
  
  #save to file
  write.table(both, file = "combined_sequences_taxa.txt", sep = "\t", quote = FALSE)
  # check that sizes of input and output are the same
  # does inputtable = newtable
  print("Checking if dimensions are the same")
  print("inputtable Rows")
  if(dim(inputtable)[1] == dim(newtable2)[1]){print("good")}
  print("inputtable Columns")
  if(dim(inputtable)[2] == dim(newtable2)[2]){print("good")}
  # does metadata = newmap
  print("metadatarows")
  if(dim(metadata)[1] == dim(newmap)[1]){print("good")}
  print("metadata columns")
  if(dim(metadata)[2] == dim(newmap)[2]){print("good")}
  # does taxa + inputtable = combined_taxa
  
  return(list(newtable = newtable, newmap = newmap, combined_taxa = both))
}
