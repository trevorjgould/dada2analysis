#' rarefaction diversity
#'
#' This function runs rarefaction stats
#'
#' @param newmap (Required).
#'  table of metadata.
#'
#' @param newtable (Required).
#' @export
#' @importFrom stats prcomp
#'  table of data.


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
scount$col <- cut(scount$scount, breaks = c(-Inf,1000,Inf), labels = c("<1000",">1000"))
p <- ggplot(scount, aes(x = scount)) + geom_dotplot(binwidth = 1000, aes(fill = col)) + geom_vline(xintercept = 1000, color = "red") + ylab("Samples") + xlab("Sequences") + scale_fill_manual(values = c("red","black"))

if (sampler < 1000){
# set minimum sequences per sample
keepers <- rowSums(newtable)>1000
newtable <- newtable[keepers,]
sampler = 1000
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

# distance:
d2.mes <- vegdist(FULL, method = "bray")
PCOA <- pcoa(d2.mes)
# from diversity.R get brayWmeta
brayWmeta$PC3pcoa <- PCOA$vectors[,3]
brayWmeta$PC2pcoa <- PCOA$vectors[,2]
brayWmeta$PC1pcoa <- PCOA$vectors[,1]

# plotting
# cage and Group NOT labelled
color4 = c("#BB3717","#54969A","#B47F6D","#0F5372")
PC1PC2 <- ggplot(brayWmeta, aes(PC1pcoa,PC2pcoa, color = as.factor(Cage))) + geom_point(size=2, aes(shape = Group)) + theme_bw() + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + scale_color_manual(values = color4,name = "Cage") + ylab(paste0("PC2: ",(var_explained[2]), "% variance")) + stat_ellipse(aes(group=as.factor(Cage), type = "norm"))
PC1PC3 <- ggplot(brayWmeta, aes(PC1pcoa,PC3pcoa, color = as.factor(Cage))) + geom_point(size=2, aes(shape = Group)) + theme_bw() + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + scale_color_manual(values = color4,name = "Cage") + ylab(paste0("PC3: ",(var_explained[3]), "% variance")) + stat_ellipse(aes(group=as.factor(Cage), type = "norm"))
#combinded_plot2 <- PC1PC2 / PC1PC3
combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
plottitle <- paste0("RAREFIED_PCoA_PC12_PC13_continuous_Group_Cage_facet_timepoint_NOT_labelled.png")
ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 12, width = 6, units = "in") 

