# load libraries
library(dplyr)

# consensus taxonomic output
# 1 given rdp and silva taxonomies from same dada2 output

# create dataframe of true/false 
tf = data.frame(matrix(nrow = nrow(silva), ncol = ncol(silva))) 
# for each column 1:7
for (x in 1:7){
tf[,x] <- silva[,x] == rdp[,x]
}
tf[is.na(tf)] <- "FALSE"

# TRUE here is consensus
consensus = data.frame(matrix(nrow = nrow(silva), ncol = ncol(silva))) 

# this is the fill in the taxa table code
# replacing NA with Unknown higher level taxa
silva$Phylum <- ifelse(is.na(silva$Phylum), paste0("Unknown ",silva$Kingdom), silva$Phylum)
silva$Class <- ifelse(is.na(silva$Class), paste0("Unknown ",silva$Phylum), silva$Class)
silva$Order <- ifelse(is.na(silva$Order), paste0("Unknown ",silva$Class), silva$Order)
silva$Family <- ifelse(is.na(silva$Family), paste0("Unknown ",silva$Order), silva$Family)
silva$Genus <- ifelse(is.na(silva$Genus), paste0("Unknown ",silva$Family), silva$Genus)
silva$Species <- ifelse(is.na(silva$Species), paste0("Unknown ",silva$Genus), silva$Species)
silva <- data.frame(lapply(silva, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
silva <- data.frame(lapply(silva, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
silva <- data.frame(lapply(silva, function(x) {gsub("Unknown Unknown", "Unknown", x)}))

# this is the fill in the taxa table code
# replacing NA with Unknown higher level taxa
rdp$Phylum <- ifelse(is.na(rdp$Phylum), paste0("Unknown ",rdp$Kingdom), rdp$Phylum)
rdp$Class <- ifelse(is.na(rdp$Class), paste0("Unknown ",rdp$Phylum), rdp$Class)
rdp$Order <- ifelse(is.na(rdp$Order), paste0("Unknown ",rdp$Class), rdp$Order)
rdp$Family <- ifelse(is.na(rdp$Family), paste0("Unknown ",rdp$Order), rdp$Family)
rdp$Genus <- ifelse(is.na(rdp$Genus), paste0("Unknown ",rdp$Family), rdp$Genus)
rdp$Species <- ifelse(is.na(rdp$Species), paste0("Unknown ",rdp$Genus), rdp$Species)
rdp <- data.frame(lapply(rdp, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
rdp <- data.frame(lapply(rdp, function(x) {gsub("Unknown Unknown", "Unknown", x)}))
rdp <- data.frame(lapply(rdp, function(x) {gsub("Unknown Unknown", "Unknown", x)}))


# if true add from cell from same position in silva 
# for each column 1:6
for (y in 1:ncol(consensus)){
  for (x in 1:nrow(consensus))
    if(tf[x,y] == "TRUE"){
      consensus[x,y] <- silva[x,y]
    }
}
colnames(consensus) <- colnames(silva)
row.names(consensus) <- row.names(silva)

# if disagreement but Not NA for both. 
# for each column 1:6
for (y in 1:ncol(consensus)){
  for (x in 1:nrow(consensus))
    if(tf[x,y] == "FALSE" && (!is.na(silva[x,y]) || !is.na(rdp[x,y]))){
      consensus[x,y] <- paste0(silva[x,y],"/",rdp[x,y])
    }
}
