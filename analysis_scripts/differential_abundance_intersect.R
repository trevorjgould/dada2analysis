# differential abundance in microbiome datasets
# no one tool has been found to accurately detect d.a.
# therefore we will do many and take the set of features at the intersection of all of them. 
# ALDEx2, ANCOM-II, corncob, MaAsLin2, metagenomeseq, DESeq2
# tools chosen based on https://www.nature.com/articles/s41467-022-28034-z
library(ALDEx2)
library(ANCOMBC)
library(Maaslin2)
library(corncob)
library(metagenomeSeq)
library(DESeq2)
library(phyloseq)
library(dplyr)
library(tidyverse)
library(purrr)
library(superheat)
library(UpSetR)
# usage:
# Rscript differential_abundance_intersect.R asvtable metadata taxonomy outputname
# assumes ASVtable is samplesAsRows
# assumes your groupings vector contains ONLY TWO GROUPS
# read in dataset
dai <- function(ASV_table, groupings, taxa, outputname){
#### testing ####
# setwd("~/Documents/xrevelo3_output/processing/")
#ASV_table = readRDS("seqtab_nochim.rds")
# groupings <- read.delim("metadata_2column.txt", row.names = 1)
#taxa <- readRDS("taxIDsilva.rds")
#################
# if standalone script
#ASV_table <- read.table(args[1], sep="\t", header=T, row.names = 1, comment.char = "", quote="", check.names = F)
# get pairwise comparisons
#groupings <- read.table(args[2], sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
# get taxa table
#taxatable = read.table(args[3], sep="\t", header=T, row.names = 1, comment.char = "", quote="", check.names = F)
# output directory
# if arg[4] exists: 
# outputname = as.character(args[4])
# else: outputname <- "output_directory"

# single column of grouping variable only
#groupings = groupings[,1]
colnames(groupings) <- "group"

# create input for x tool
print("creating tables")
samplesAsRows <- ASV_table
samplesAsColumns <- as.data.frame(t(ASV_table))
phyloseqObject <- phyloseq(otu_table(ASV_table, taxa_are_rows = FALSE), sample_data(groupings), tax_table(taxa))
metagenomeSeqObject <- phyloseq_to_metagenomeSeq(phyloseqObject)
ancombcobject = mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseqObject)
# run ALDEx2
  # https://www.bioconductor.org/packages/devel/bioc/manuals/ALDEx2/man/ALDEx2.pdf
  # reads: A non-negative, integer-only data.frame or matrix with unique names for all rows and columns. Rows should contain genes and columns should contain sequencing read counts (i.e., sample vectors). Rows with 0 reads in each sample are deleted prior to analysis.
  # conditions: A character vector. A description of the data structure used for testing. Typically, a vector of group labels. For aldex.glm, use a model.matrix.
print("running ALDEx2")
ALDEx2_results <- aldex(reads=samplesAsColumns, conditions = groupings$group, mc.samples = 128, 
                        test="t", effect=TRUE, include.sample.summary = FALSE, 
                        verbose=T, denom="all")
filename1 <- paste0(outputname,"_ALDEx2_result.txt")
write.table(ALDEx2_results, file=filename1, quote=FALSE, sep='\t', col.names = NA)
print("finished ALDEx2")
# run ANCOMBC
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
# phyloseq input
# currently running ancombc2
print("running ANCOMBC")
ANCOMBC_results <- ancombc2(data = ancombcobject, struc_zero = "TRUE", group = colnames(groupings[1]), 
                            alpha=0.05, p_adj_method = "holm", fix_formula = colnames(groupings[1]))
filename2 <- paste0(outputname,"_ANCOMBC_result.txt")
write.table(ANCOMBC_results$res, file=filename2, quote=FALSE, sep='\t', col.names = NA)
print("finished ANCOMBC")
# run Maaslin2
# input_data: tab delimited, features as columns, samples as rows
# input_metadata: tab delimited, metadata as columns, samples as rows
print("running Maaslin2")
Maaslin2_results <- Maaslin2(samplesAsRows, groupings, outputname, transform = "LOG", 
                             fixed_effects = c(colnames(groupings[1])), standardize = FALSE, 
                             plot_heatmap = F, plot_scatter = F)

filename3 <- paste0(outputname,"_Maaslin2_result.txt")
write.table(Maaslin2_results$results, file=filename3, quote=FALSE, sep='\t', col.names = NA)
print("finished Maaslin2")
# run corncob
# phyloseq input
# https://cran.r-project.org/web/packages/corncob/corncob.pdf
# https://rdrr.io/github/bryandmartin/corncob/man/differentialTest.html
print("running corncob")
my_formula <- as.formula(paste("~",colnames(groupings[1]),sep=" ", collapse = ""))
corncob_results <- corncob::differentialTest(formula= my_formula,
                                     phi.formula = my_formula, phi.formula_null = my_formula,
                                     formula_null = ~ 1, test="Wald", data=phyloseqObject,
                                     fdr_cutoff = 0.05)
filename4 <- paste0(outputname,"_corncob_result.txt")

ccfdr <- as.data.frame(corncob_results$p_fdr)
ccfdr$p = corncob_results$p
colnames(ccfdr) <- c("Pvaluefdr","Pvalue")
write.table(ccfdr, file=filename4, quote=FALSE, sep='\t', col.names = NA)
print("finished corncob")
# run metagenomeSeq
# data: rows = features, columns = samples
# takes metagenomeSeqObject
print("running metagenomeSeq")
p <- cumNormStatFast(metagenomeSeqObject)
metagenomeSeqObject2 <- cumNorm(metagenomeSeqObject, p = p)
pd <- pData(metagenomeSeqObject2)
mod <- model.matrix( ~ 1 + group, data = pd)
regres <- fitFeatureModel(metagenomeSeqObject2, mod)
metagenomeSeq_result <- MRfulltable(regres, number = nrow(featureData(metagenomeSeqObject)))
filename5 <- paste0(outputname,"_metagenomeSeq_result.txt")
write.table(metagenomeSeq_result, file=filename5, quote=FALSE, sep='\t', col.names = NA)
print("finished metagenomeSeq")
# run DESeq2
# columns: samples
print("running deseq2")
groupings <- groupings[order(row.names(groupings)), , drop = FALSE]
dds <- DESeq2::DESeqDataSetFromMatrix(countData = samplesAsColumns, colData=groupings, design =~group)
dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")
res <- results(dds_res, tidy=T, format="DataFrame")
rownames(res) <- res$row
deseq2_result <- res[,-1]
filename6 <- paste0(outputname,"_deseq2_result.txt")
write.table(deseq2_result, file=filename6, quote=FALSE, sep="\t", col.names = NA)
print("finished deseq2")
# We now have 6 tables of results that need to be combined
# ALDEx2_result
# ANCOMBC_result
# Maaslin2_result
# corncob_result
# metagenomeSeq_result
# deseq2_result
print("combining methods")
results_tables <- list.files(path = ".", pattern = "*_result.txt")
if (length(results_tables)<6){
  message1 = paste0("Only found ",length(results_tables)," results tables. Some methods are missing results.")
  print(message1)
} else {
  message1 = paste0("Found ",length(results_tables)," results tables. proceeding to merge.")
  print(message1)}

# for each table pull out pval and/or adjusted pval / effect. 
# ALDEx2_result
df1 <- as.data.frame(cbind(row.names(ALDEx2_results),ALDEx2_results$effect))# greater than 1 is recommended as reproducible significance cutoff
colnames(df1) <- c("seqid","ALDEx2_adj")
# ANCOMBC_result
df2 <- as.data.frame(cbind(ANCOMBC_results$res$taxon,dplyr::select(ANCOMBC_results$res, starts_with("p_"), -starts_with("p_(")),dplyr::select(ANCOMBC_results$res, starts_with("q_"), -starts_with("q_("))))
colnames(df2) <- c("seqid","ANCOMBC2_Pvalue","ANCOMBC2_adj")
# Maaslin2_result
df3 <- as.data.frame(cbind(Maaslin2_results$results$feature,Maaslin2_results$results$pval,Maaslin2_results$results$qval))
colnames(df3) <- c("seqid","Maaslin2_Pvalue","Maaslin2_adj")
# corncob_result
df4 <- as.data.frame(cbind(row.names(ccfdr),ccfdr$Pvalue,ccfdr$Pvaluefdr))
colnames(df4) <- c("seqid","corncob_Pvalue","corncob_adj")
# metagenomeSeq_result
df5 <- as.data.frame(cbind(row.names(metagenomeSeq_result),metagenomeSeq_result$pvalues,metagenomeSeq_result$adjPvalues))
colnames(df5) <- c("seqid","metagenomeseq_Pvalue","metagenomeseq_adj")
# deseq2_result
df6 <- as.data.frame(cbind(row.names(deseq2_result),deseq2_result$pval,deseq2_result$padj))
colnames(df6) <- c("seqid","desq2_Pvalue","deseq2_adj")
# Merge 6 tables together
dflist <- list(df1, df2, df3, df4, df5, df6)
dfall <- dflist %>% purrr::reduce(full_join, by='seqid')
filename7 <- paste0(outputname,"_merged_all_result.txt")
write.table(dfall, file=filename7, quote=FALSE, sep='\t', col.names = NA)
print("calculating results")
# next we will take the combined result of that function "dfall" 
# and present results 
# split full table by Pvalue and adj Pvalue
adjustedOnly <- dplyr::select(dfall, ends_with("_adj"))
PvalueOnly <- dplyr::select(dfall, ends_with("_Pvalue"))
row.names(adjustedOnly) <- dfall$seqid
row.names(PvalueOnly) <- dfall$seqid
write.table(PvalueOnly, file = "PvalueOnly.txt", sep = "\t", quote = FALSE)
write.table(adjustedOnly, file = "adjustedOnly.txt", sep = "\t", quote = FALSE)
adjustedOnly <- read.table("adjustedOnly.txt", header = TRUE, row.names = 1)
PvalueOnly <- read.table("PvalueOnly.txt", header = TRUE, row.names = 1)

# heatmaps but not terribly useful in larger datasets
print("making heatmaps")
p2 <- filter(PvalueOnly, rowSums(is.na(PvalueOnly))<=1)
a2 <- filter(adjustedOnly, rowSums(is.na(adjustedOnly))<=1)
superheat(p2, left.label = "none", pretty.order.rows = TRUE, pretty.order.cols = TRUE, heat.pal = c("#CC0000","#FFCC66","#336699"),heat.pal.values = c(0,0.05,1))
superheat(a2[,2:6], left.label = "none", pretty.order.rows = TRUE, pretty.order.cols = TRUE, heat.pal = c("#CC0000","#FFCC66","#336699"),heat.pal.values = c(0,0.05,1))
# just for any tool returning p<0.05
p2 <- filter_all(PvalueOnly, any_vars(. < 0.05))
a2 <- filter_all(adjustedOnly[,2:6], any_vars(. < 0.05))
write.table(p2, file = "PvalueOnly_Any_tool_significant_only.txt", sep = "\t", quote = FALSE)
write.table(a2, file = "adjustedOnly_Any_tool_significant_only.txt", sep = "\t", quote = FALSE)
p2 <- filter_all(PvalueOnly, any_vars(. < 0.05))
superheat(p2, left.label = "none", pretty.order.rows = TRUE, pretty.order.cols = TRUE, heat.pal = c("#CC0000","#FFCC66","#336699"),heat.pal.values = c(0,0.05,1))
superheat(a2, left.label = "none", pretty.order.rows = TRUE, pretty.order.cols = TRUE, heat.pal = c("#CC0000","#FFCC66","#336699"),heat.pal.values = c(0,0.05,1))

# get only p value <= 0.05 more than 1 tool
PvalueOnly$count <- rowSums(PvalueOnly <= 0.05)
ConsensusPvalue <- subset(PvalueOnly, count > 1)
ConsensusPvalue <- dplyr::select(ConsensusPvalue, -ends_with("count"))
write.table(ConsensusPvalue, file = "ConsensusSigPvalueOnly.txt", sep = "\t", quote = FALSE)

adjustedOnly$count <- rowSums(adjustedOnly <= 0.05)
ConsensusAdjustedOnly <- subset(adjustedOnly, count > 1)
ConsensusAdjustedOnly <- dplyr::select(ConsensusAdjustedOnly, -ends_with("count"))
write.table(ConsensusAdjustedOnly, file = "ConsensusSigadjustedOnly.txt", sep = "\t", quote = FALSE)
vectorOfTables <- vector(mode = "list", length = 2)
vectorOfTables[[1]] <- ConsensusPvalue
vectorOfTables[[2]] <- ConsensusAdjustedOnly
names(vectorOfTables) <- c("ConsensusPvalue","ConsensusAdjustedOnly")
print("finished")
return(vectorOfTables)
}

# example
vot <- dai(t1b, t2b, t3b, outputdir)
# From there we could take the significant ASVs and plot the counts per sample vs Group
both <- merge(asv2,groupings, by = "row.names")
both <- column_to_rownames(both, var = "Row.names")
melted <- melt(both)
ggplot(melted, aes(Group,log(value))) + geom_boxplot() + facet_wrap(~variable) + theme_bw() + ylab("log(count)")

library(UpSetR)
# create UpSetR format table
# here we need a table of 1 or 0 for sig/not 
getUpset <- function(intab,cutoff){
  oldrn <- rownames(intab)
  row.names(intab) <- 1:nrow(intab)
  intab2 <- intab
  intab2[intab<=cutoff] = 1
  intab2[intab>cutoff] = 0
  intab2[is.na(intab2)] = 0
  pdf(file="upsetR_plot.pdf",onefile=FALSE) # or other device
  plot1 <- upset(intab2, order.by = "freq")
  print(plot1)
  dev.off()
  return(intab2)
}
gu <- getUpset(vectorOfTables[[1]],0.05)

# same thing for the adjusted table with one addition
getUpsetAdjust <- function(intab,cutoff,aldecut){
  oldrn <- rownames(intab)
  row.names(intab) <- 1:nrow(intab)
  intab2 <- intab
  intab2[intab<=cutoff] = 1
  intab2[intab>cutoff] = 0
  intab2[is.na(intab2)] = 0
  intab2[,1] = intab[,1]
# aldex2 is not a p-value. greater than 1 is likely to be reproducible effect
# here we set make the cutoff flexible 
  good = intab2$ALDEx2_adj > aldecut
  if(any(good == TRUE) > 0){
  good["TRUE"] = 1
  good["FALSE"] = 0
  intab2$ALDEx2_adj = good[1:6]
  }
  pdf(file="Adjusted_upsetR_plot.pdf",onefile=FALSE) # or other device
  intab3 <- intab2[, which(colSums(intab2) != 0)]
  plot2 <- upset(intab3, order.by = "freq", nsets = 6)
  print(plot2)
  dev.off()
  rownames(intab2)<- oldrn
  return(intab2)
}
gua <- getUpsetAdjust(vectorOfTables[[2]],0.05,1)
