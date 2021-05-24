#combined functions
newcolors <- c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e','#003c30')
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
#save to file
write.table(both, file = "combined_sequences_taxa.txt", sep = "\t", quote = FALSE)
return(list(newtable = newtable, newmap = newmap, combined_taxa = both))
}

# read in files
Diversity_Plots <- function(brayWmeta,newmap){
# brayWmeta <- read.table('proportional_diversity_stats.txt')
# newmap <- read.table("Metadata_common.txt")
var_explained = (brayWmeta$EV/sum(brayWmeta$EV))*100
var_explained = format(round(var_explained, 2), nsmall = 2)

Adiv <- function(x) {
    ShanD <-  ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$shannon, fill = x)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white")+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + scale_discrete_manual(values = newcolors, aesthetics = "fill")
	SimD <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$simpson, fill = x)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white")+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Simpson") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + scale_discrete_manual(values = newcolors, aesthetics = "fill")
	SimI <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$shannon, fill = x)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white")+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("invsimpson") + scale_discrete_manual(values = newcolors, aesthetics = "fill")
    #combinded_plot <- ShanD / SimD / SimI
    combinded_plot <- gridExtra::grid.arrange(ShanD,SimD,SimI, ncol = 1)
    plottitle <- paste0("AlphaDiversity_plots_",x,".png")
    ggplot2::ggsave(combinded_plot, file=plottitle, dpi=800, height = 12, width = 6, units = "in")
}

# ggplot functions for PCoAs
bdiv <- function(j) {
    PC1PC2 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC2, colour = j)) + ggplot2::geom_point(pch=21, size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC2: ",(var_explained[2]), "% variance")) + scale_color_manual(values=newcolors)
    PC1PC3 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC3, colour = j)) + ggplot2::geom_point(pch=21, size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC3: ",(var_explained[3]), "% variance")) + scale_color_manual(values=newcolors)
    #combinded_plot2 <- PC1PC2 / PC1PC3
    combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
    plottitle <- paste0("PCoA_PC12_PC13_continuous",j,".png")
    ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
}

    # make plots
    num_cols <- unlist(lapply(outtab$newmap, is.numeric))   
    fac_cols <- !num_cols
    lapply(colnames(outtab$newmap), Adiv)
    lapply(colnames(outtab$newmap), bdiv)
}

# metadata table
diversity <- function(newmap,newtable){
#newmap <- read.table("Metadata_common.txt")
#newtable <- read.table("Sequence_table_common.rds")

# get counts
newmap$count <- rowSums(newtable)

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

#PCA
d.mes <- prcomp(data.CLR, scale = FALSE)
var_explained = (d.mes$sdev^2/sum(d.mes$sdev^2))*100
var_explained = format(round(var_explained, 2), nsmall = 2)

#eigenvalues
newmap$EV <- d.mes$sdev^2
brayWmeta <- cbind(newmap,d.mes$x[,1:4])

#alpha_diversity_stats
propdist <- sweep(newtable, 1, rowSums(newtable),'/')
brayWmeta$shannon <- vegan::diversity(propdist, index = "shannon")
brayWmeta$simpson <- vegan::diversity(propdist, index = "simpson")
brayWmeta$invsimpson <- vegan::diversity(propdist, index = "invsimpson")
write.table(brayWmeta, file="proportional_diversity_stats.txt", quote = FALSE)
return(brayWmeta)
}

# reads in table from Make_Tables.R
Make_Taxa_Tables <- function(x){
combined_taxa <- read.table(file = x, sep = "\t")
# taxonomy_tables
# files to use for taxa:
# metadata table
#metadata <- read.table("Metadata_common.txt")
# taxa
#newtable <- read.table("combined_sequences_taxa.txt", sep = "\t", check.names = FALSE)

#make split taxa tables
levels <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
n <- (ncol(combined_taxa) - 7)
`%>%` <- dplyr::`%>%`
Domain_table <- combined_taxa %>% dplyr::select(Kingdom,1:n)
Phylum_table <- combined_taxa %>% dplyr::select(Phylum,1:n)
Class_table <- combined_taxa %>% dplyr::select(Class,1:n)
Order_table <- combined_taxa %>% dplyr::select(Order,1:n)
Family_table <- combined_taxa %>% dplyr::select(Family,1:n)
Genus_table <- combined_taxa %>% dplyr::select(Genus,1:n)
Species_table <- combined_taxa %>% dplyr::select(Species,1:n)
KT <- plyr::ddply(Domain_table, "combined_taxa$Kingdom", plyr::numcolwise(sum))
PT <- plyr::ddply(Phylum_table, "combined_taxa$Phylum", plyr::numcolwise(sum))
CT <- plyr::ddply(Class_table, "combined_taxa$Class", plyr::numcolwise(sum))
OT <- plyr::ddply(Order_table, "combined_taxa$Order", plyr::numcolwise(sum))
FT <- plyr::ddply(Family_table, "combined_taxa$Family", plyr::numcolwise(sum))
GT <- plyr::ddply(Genus_table, "combined_taxa$Genus", plyr::numcolwise(sum))
ST <- plyr::ddply(Species_table, "combined_taxa$Species", plyr::numcolwise(sum))
KT = setNames(data.frame(t(KT[,-1])), KT[,1])
PT = setNames(data.frame(t(PT[,-1])), PT[,1])
CT = setNames(data.frame(t(CT[,-1])), CT[,1])
OT = setNames(data.frame(t(OT[,-1])), OT[,1])
FT = setNames(data.frame(t(FT[,-1])), FT[,1])
GT = setNames(data.frame(t(GT[,-1])), GT[,1])
ST = setNames(data.frame(t(ST[,-1])), ST[,1])
write.table(KT, file = "Kingdom_taxonomy.txt", quote = FALSE, sep = "\t")
write.table(PT, file = "Phylum_taxonomy.txt", quote = FALSE, sep = "\t")
write.table(CT, file = "Class_taxonomy.txt", quote = FALSE, sep = "\t")
write.table(OT, file = "Order_taxonomy.txt", quote = FALSE, sep = "\t")
write.table(FT, file = "Family_taxonomy.txt", quote = FALSE, sep = "\t")
write.table(GT, file = "Genus_taxonomy.txt", quote = FALSE, sep = "\t")
write.table(ST, file = "Species_taxonomy.txt", quote = FALSE)
return(list(KT=KT,PT=PT,CT=CT,OT=OT,FT=FT,GT=GT,ST=ST))
}

Taxonomy_Plots <- function(meta){
# read in tables
#meta <- read.table("Metadata_common.txt", sep = "\t", check.names = FALSE)
KT <- read.table(file = "Kingdom_taxonomy.txt", sep = "\t", check.names = FALSE)
PT <- read.table(file = "Phylum_taxonomy.txt", sep = "\t", check.names = FALSE)
CT <- read.table(file = "Class_taxonomy.txt", sep = "\t", check.names = FALSE)
OT <- read.table(file = "Order_taxonomy.txt", sep = "\t", check.names = FALSE)
FT <- read.table(file = "Family_taxonomy.txt", sep = "\t", check.names = FALSE)
GT <- read.table(file = "Genus_taxonomy.txt", sep = "\t", check.names = FALSE)
ST <- read.table(file = "Species_taxonomy.txt", sep = "\t", check.names = FALSE)

# Here we are taking the taxonomy tables created above
# taxa with less than a sum of 0.1 total proportion over all samples are merged
# into a column of OTHER with the NA column if it is present.
# NEXT IS TO DYNAMICALLY SET THE 0.1 TO LEAVE 20 TAXA TOTAL WITH OTHER AS REST.
taxa_list <- list(KT,PT,CT,OT,FT,GT,ST)
taxa_names <- list("Kingdom","Phylum","Class","Order","Family","Genus","Species")
for (x in 2:7){
XT_other = 0
xt=0
XTtax=0
name_label <- taxa_names[x]
xt=as.data.frame(taxa_list[x])
xt=xt[(rowSums(xt)>0),]

# fix rownames . to -
row.names(meta) <- gsub("-", "\\.", row.names(meta))

#find common names
common <- intersect(rownames(meta),rownames(xt))
# get just the overlapping samples
meta <- meta[common,, drop = FALSE]
xt <- xt[common,, drop = FALSE]


XTtax = sweep(xt, 1, rowSums(xt),'/')
# KEEP TOP TAXA AND SORT REST INTO OTHER CATEGORY AND PLOT
# get taxa that are present > 0.1
XT_other <- XTtax[,(colSums(XTtax)>=0.1)]
# sum taxa less then 0.1

if (table(colSums(XTtax)<0.1)["TRUE"] > 1){
XT_other$other <- rowSums(XTtax[,(colSums(XTtax)<0.1)])

#
if("NA." %in% colnames(XT_other))
{
cat("Merging NA column with Other!\n");
XT_other$Other <- XT_other$'NA.' + XT_other$other
# remove columns NA and other
XT_other <- XT_other[ , -which(names(XT_other) %in% c("NA.","other"))]
}
b <- sum(nrow(XT_other))
a <- as.integer(sum(rowSums(XT_other)))

if(a==b){
cat("taxa are looking good\n")
}
if(a!=b){
 cat("taxa do not sum to 100. there is something wrong\n")
}
}


both <- cbind(XT_other,meta)
both$Samples <- row.names(meta)
melted <- reshape2::melt(both, id.vars = c(colnames(meta),"Samples"))
filename <- paste0(name_label,"_taxonomy_other.png")
filename2 <- paste0(name_label,"_taxonomy_other.txt")
p <- ggplot2::ggplot(melted, ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + ggplot2::scale_fill_discrete(name = name_label)
ggplot2::ggsave(p, file = filename, dpi  = 800, width = 10, height = 8, units = "in")
write.table(both, file = filename2, sep = "\t", quote = FALSE)
#V1p <- p + ggplot2::facet_wrap(.~ var1, scales = "free_x")
#filename = paste0(name_label,"-",x,"-V1_split_taxonomy_other.png")
#ggplot2::ggsave(V1p, file = filename, dpi  = 800, width = 10, height = 8, units = "in")

}
}

# input dada2 sequence table
t1 <- readRDS('seqtab_nochim.rds')
# input metadata
t2 <- read.table('metadata.txt', sep = '\t', comment='', head=TRUE, row.names=1, check.names = FALSE)
# input taxonomy
t3 <- readRDS("taxID.rds")

outtab <- Create_Tables(t1,t2,t3)
combined_taxa <- read.table(file = "combined_sequences_taxa.txt", sep = "\t")

taxa_out <- Make_Taxa_Tables("combined_sequences_taxa.txt")
Taxonomy_Plots(outtab$newtable)
brayWmeta <- diversity(outtab$newmap,outtab$newtable)
Diversity_Plots(brayWmeta, outtab$newmap)