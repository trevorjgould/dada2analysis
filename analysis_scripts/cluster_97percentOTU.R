# module load R/4.3.0-openblas
library(tibble)
library(dplyr)
# Packages that are required but not loaded:
# library(DECIPHER)
# library(Biostrings)

nproc <- 24 # set to number of cpus/processors to use for the clustering
seqtab <- readRDS("seqtab_nochim.rds")
taxa <- readRDS("taxa.rds")

asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

## Find clusters of ASVs to form the new OTUs
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)

clusters <- DECIPHER::TreeLine(
  myDistMatrix=d,
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU
  type = "clusters",
  processors = nproc)
  
  
## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters <- clusters %>% add_column(sequence = asv_sequences)

merged_seqtab <- seqtab %>% 
  t %>%
  rowsum(clusters$cluster) %>%
  t
  
# Optional renaming of clusters to OTU<cluster #>
colnames(merged_seqtab) <- paste0("OTU", colnames(merged_seqtab))

saveRDS(merged_seqtab,file = "clustered_ASV.rsv")
write.table(merged_seqtab, file = "clustered_ASV.txt", sep = "\t", quote = FALSE, col.names=NA)

both_clust <- cbind(taxa,clusters)
write.table(both_clust, file = "clusters_taxa.txt", sep = "\t", quote = FALSE, row.names = FALSE)