#' Diversity Stats
#'
#' This function runs stats for alpha and beta diversity
#'
#' @param newmap table
#' @param newtable table
#'
#' @export

diversity <- function(newmap,newtable){
#newmap <- read.table("Metadata_common.txt")
#newtable <- read.table("Sequence_table_common.rds")

# get counts
newmap$count <- rowSums(newtable)
#newmap <- newmap[-which(newmap$Eligibility_group == ""), ]

# 0 out NA cells
newtable[is.na(newtable)] <- 0
newtable2 <- t(newtable)

data.CLR <- clr_transform(newtable2)
data.CLR <- t(data.CLR)

# get just the overlapping samples
#meta <- newmap[-which(newmap$Eligibility_group == ""), ]
common <- intersect(rownames(data.CLR),rownames(newmap))
data.CLR <- data.CLR[common,, drop = FALSE]
newmap <- newmap[common,, drop = FALSE]
newtable <- newtable[common,, drop = FALSE]

#PCA
d.mes <- stats::prcomp(data.CLR, scale = FALSE)
var_explained = (d.mes$sdev^2/sum(d.mes$sdev^2))*100
var_explained = format(round(var_explained, 2), nsmall = 2)

#eigenvalues
newmap$EV <- d.mes$sdev^2
brayWmeta <- cbind(newmap,d.mes$x[,1:4])
newtable2 <- t(newtable)
brayWmeta$chao1 <- apply(newtable2, 2, OTUtable::chao1)
#alpha_diversity_stats
propdist <- sweep(newtable, 1, rowSums(newtable),'/')
brayWmeta$shannon <- vegan::diversity(propdist, index = "shannon")
brayWmeta$simpson <- vegan::diversity(propdist, index = "simpson")
#brayWmeta$invsimpson <- vegan::diversity(propdist, index = "invsimpson")
return(brayWmeta)
}
######################################################
