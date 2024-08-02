#biplot
library(ggplot2)
library(factoextra)
library(reshape2)
library(dplyr)

# input
# taxa <- read.delim("Genus_taxonomy_other.txt")
# meta <- read.delim("metadata.txt", row.names=1, check.names = FALSE)

# makes a biplot for each column in meta based on taxa table.
biplot_function <- function(taxa,meta){
taxatable = as.data.frame(taxa)
metadata = as.data.frame(meta)
n = ncol(taxa)
n = n-7
G2 = taxatable[,1:n]
mod = prcomp(na.omit(G2), scale=TRUE)

for (x in (1:ncol(meta))){
p <- fviz_pca_biplot(mod, label="var", repel=TRUE, habillage = meta[,x], select.var = list(cos2 = 0.25),ggtheme = theme_bw())
f1 <- colname(meta[x])
f2 <- paste0("genus_biplot_",f1,".png")
ggsave(p, file = f2, dpi=800)
}
}

# biplot_function(taxa,meta)