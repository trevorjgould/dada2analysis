#' Diversity Stats
#'
#' This function runs stats for alpha and beta diversity
#'
#' @param newmap (Required).
#'  table of metadata.
#'
#' @param newtable (Required).
#' @export
#' @importFrom stats prcomp
#'  table of data.

# metadata table
diversity <- function(newmap,newtable){
newmap <- read.table("Metadata_common.txt")
newtable <- read.table("Sequence_table_common.rds")

# get counts
newmap$count <- rowSums(newtable)
#newmap <- newmap[-which(newmap$Eligibility_group == ""), ]

# 0 out NA cells
newtable[is.na(newtable)] <- 0
newtable2 <- t(newtable)
## clr + imputation function
## Borrowed from Gabe
## Note, this assumes taxa as rows and samples as columns
clr_transform <- function(x){
  clr.taxa <- x
  clr.taxa = t(clr.taxa); eps = 0.5
  clr.taxa = clr.taxa*(1 -rowSums(clr.taxa==0)*eps/rowSums(clr.taxa))
  clr.taxa[clr.taxa==0]=eps
  clr.taxa = sweep(clr.taxa,1,rowSums(clr.taxa),'/');
  ls = log(clr.taxa)
  clr.taxa = t(ls - rowMeans(ls))
  clr.taxa = clr.taxa[,!is.nan(colSums(clr.taxa))]
  return(clr.taxa)
}

data.CLR <- clr_transform(newtable2)
data.CLR <- t(data.CLR)

# get just the overlapping samples
#meta <- newmap[-which(newmap$Eligibility_group == ""), ]
common <- intersect(rownames(data.CLR),rownames(newmap))
data.CLR <- data.CLR[common,, drop = FALSE]
newmap <- newmap[common,, drop = FALSE]
newtable <- newtable[common,, drop = FALSE]

#PCA
d.mes <- prcomp(data.CLR, scale = FALSE)
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


############
# RAREFEACTION
library(vegan)
library(plyr)
library(dplyr)
# first we need to see what the sample counts look like to pick a sample number
# get counts
# This assumes samples = ROWS
min(rowSums(newtable))
hist(rowSums(newtable))
sampler = min(rowSums(newtable))
scount <- as.data.frame(rowSums(newtable))
print(paste("lowest sample size is: ", sampler))
if (sampler < 1000){
  print("Setting rarefaction threshold to: 1000")
  # set minimum sequences per sample
  keepers <- rowSums(newtable)>1000
  newtable <- newtable[keepers,]
  sampler = 1000
  scount$col <- cut(scount$scount, breaks = c(-Inf,1000,Inf), labels = c("<1000",">1000"))
  ggplot(scount, aes(x = scount)) + geom_dotplot(binwidth = 1000, aes(fill = col)) + geom_vline(xintercept = 1000, color = "red") + ylab("Samples") + xlab("Sequences") + scale_fill_manual(values = c("red","black"))
} else {
  print(paste("Setting rarefaction threshold to: ", sampler))
}

# then we rarefy the datatable at chosen sample size
rarefyAt <- function(newtable,sampler){
  # we want to do rrarefy 1000 times and get the average   
  # lets try with 10 first
  res <- lapply(as.list(1:1000), function(x) rrarefy(newtable, sample=sampler))
  FULL <- aaply(laply(res, as.matrix), c(2, 3), mean)
  rm(res)
  return(FULL)
}
out <- rarefyAt(newtable,sampler)
# now we're going to compare the 1000 times rarefied table to the original run via clr transform
propdistout <- sweep(out, 1, rowSums(out),'/')
# distance:
d2.mes <- vegdist(propdistout, method = "bray")
PCOA <- cmdscale(d2.mes, k = 3, eig = TRUE)
# from diversity.R get brayWmeta
brayWmeta$PC3pcoa <- PCOA$points[,3]
brayWmeta$PC2pcoa <- PCOA$points[,2]
brayWmeta$PC1pcoa <- PCOA$points[,1]
###################
round(PCOA$eig*100/sum(PCOA$eig),2)
write.table(brayWmeta, file="proportional_diversity_stats.txt", quote = FALSE)
outlist <- list("brayWmeta"=brayWmeta,"PCOA" = PCOA)
return(outlist)
}

######################################################
adonis_stats <- function(newmap,newtable,var1,var2){
  set.seed(12345)
  # get counts
  newmap$Scount <- rowSums(newtable)
  
  # 0 out NA cells
  newtable[is.na(newtable)] <- 0
  newtable2 <- t(newtable)
  
  ## clr + imputation function
  ## Borrowed from Gabe
  ## Note, this assumes taxa as rows and samples as columns
  clr_transform <- function(x){
    clr.taxa <- x
    clr.taxa = t(clr.taxa); eps = 0.5
    clr.taxa = clr.taxa*(1 -rowSums(clr.taxa==0)*eps/rowSums(clr.taxa))
    clr.taxa[clr.taxa==0]=eps
    clr.taxa = sweep(clr.taxa,1,rowSums(clr.taxa),'/');
    ls = log(clr.taxa)
    clr.taxa = t(ls - rowMeans(ls))
    clr.taxa = clr.taxa[,!is.nan(colSums(clr.taxa))]
    return(clr.taxa)
  }
  data.CLR <- clr_transform(newtable2)
  data.CLR <- t(data.CLR)
  
  # get just the overlapping samples
  common <- intersect(rownames(data.CLR),rownames(newmap))
  data.CLR <- data.CLR[common,, drop = FALSE]
  newmap <- newmap[common,, drop = FALSE]
  newtable <- newtable[common,, drop = FALSE]
  
  # dist
  out <- dist(data.CLR)
  
  # adonis
  aout <- adonis2(out ~ var1 + Scount, data = newmap, permutations = 10000)
  return(aout)
}
adonis_stats(newmap,newtable,var1,var2)