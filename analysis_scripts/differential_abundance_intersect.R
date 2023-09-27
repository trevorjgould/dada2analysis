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
# usage:
# Rscript differential_abundance_intersect.R asvtable metadata taxonomy outputname
# assumes ASVtable is samplesAsRows
# read in dataset
dai <- function(ASV_table, groupings, taxa, outputname){
#### testing ####
# setwd("~/Documents/xrevelo3_output/processing/")
# ASV_table = readRDS("seqtab_nochim.rds")
# groupings <- read.delim("metadata_2column.txt", row.names = 1)
# taxa <- readRDS("taxIDsilva.rds")
# studyname <- "output_dir"
#################
# if standalone script
#ASV_table <- read.table(args[1], sep="\t", header=T, row.names = 1, comment.char = "", quote="", check.names = F)
# get pairwise comparisons
#groupings <- read.table(args[2], sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
# get taxa table
#taxatable = read.table(args[3], sep="\t", header=T, row.names = 1, comment.char = "", quote="", check.names = F)
# output directory
# if arg[4] exists: 
# studyname = as.character(args[4])
# else: studyname <- "output_directory"

# single column of grouping variable only
groupings = groupings[,1]

# create input for x tool
samplesAsRows <- ASV_table
samplesAsColumns <- as.data.frame(t(ASV_table))
phyloseqObject <- phyloseq(otu_table(ASV_table, taxa_are_rows = FALSE), sample_data(groupings), tax_table(taxa))
metagenomeSeqObject <- phyloseq_to_metagenomeSeq(phyloseqObject)
ancombcobject = mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseqObject)
# run ALDEx2
  # https://www.bioconductor.org/packages/devel/bioc/manuals/ALDEx2/man/ALDEx2.pdf
  # reads: A non-negative, integer-only data.frame or matrix with unique names for all rows and columns. Rows should contain genes and columns should contain sequencing read counts (i.e., sample vectors). Rows with 0 reads in each sample are deleted prior to analysis.
  # conditions: A character vector. A description of the data structure used for testing. Typically, a vector of group labels. For aldex.glm, use a model.matrix.
ALDEx2_results <- aldex(reads=samplesAsColumns, conditions = groupings[,1], mc.samples = 128, 
                        test="t", effect=TRUE, include.sample.summary = FALSE, 
                        verbose=T, denom="all")
filename1 <- paste0(studyname,"_ALDEx2_result.txt")
write.table(ALDEx2_results, file=filename1, quote=FALSE, sep='\t', col.names = NA)

# run ANCOMBC
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
# phyloseq input
# currently running ancombc2
ANCOMBC_results <- ancombc2(data = ancombcobject, struc_zero = "TRUE", group = colnames(groupings[1]), 
                            alpha=0.05, p_adj_method = "holm", fix_formula = colnames(groupings[1]))
filename2 <- paste0(studyname,"_ANCOMBC_result.txt")
write.table(ANCOMBC_results$res, file=filename2, quote=FALSE, sep='\t', col.names = NA)

# run Maaslin2
# input_data: tab delimited, features as columns, samples as rows
# input_metadata: tab delimited, metadata as columns, samples as rows
Maaslin2_results <- Maaslin2(samplesAsRows, groupings, outputname, transform = "LOG", 
                             fixed_effects = c(colnames(groupings[1])), standardize = FALSE, 
                             plot_heatmap = F, plot_scatter = F)

filename3 <- paste0(studyname,"_Maaslin2_result.txt")
write.table(Maaslin2_results$results, file=filename3, quote=FALSE, sep='\t', col.names = NA)

# run corncob
# phyloseq input
# https://cran.r-project.org/web/packages/corncob/corncob.pdf
# https://rdrr.io/github/bryandmartin/corncob/man/differentialTest.html
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

# run metagenomeSeq
# data: rows = features, columns = samples
# takes metagenomeSeqObject
p <- cumNormStatFast(metagenomeSeqObject)
metagenomeSeqObject2 <- cumNorm(metagenomeSeqObject, p = p)
pd <- pData(metagenomeSeqObject2)
mod <- model.matrix( ~ 1 + colnames(groupings[1]), data = pd)
regres <- fitFeatureModel(metagenomeSeqObject2, mod)
metagenomeSeq_result <- MRfulltable(regres, number = nrow(featureData(metagenomeSeqObject)))
filename5 <- paste0(outputname,"_metagenomeSeq_result.txt")
write.table(metagenomeSeq_result, file=filename5, quote=FALSE, sep='\t', col.names = NA)

# run DESeq2
# columns: samples
dds <- DESeq2::DESeqDataSetFromMatrix(countData = samplesAsColumns,
                                      colData=groupings,
                                      design = ~ colnames(groupings[1]))
dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")
res <- results(dds_res, tidy=T, format="DataFrame")
rownames(res) <- res$row
deseq2_result <- res[,-1]
filename6 <- paste0(outputname,"_deseq2_result.txt")
write.table(deseq2_result, file=filename6, quote=FALSE, sep="\t", col.names = NA)

# We now have 6 tables of results that need to be combined
# ALDEx2_result
# ANCOMBC_result
# Maaslin2_result
# corncob_result
# metagenomeSeq_result
# deseq2_result

results_tables <- list.files(path = ".", pattern = "*_result.txt")
if (length(results_tables)<6){
  print0("Only found ",length(results_tables)," results tables. Something went wrong.")
} else {print0("Found ",length(results_tables)," results tables. proceeding to merge.")}

# for each table pull out pval and/or adjusted pval / effect. 
# ALDEx2_result
df1 <- as.data.frame(cbind(row.names(ALDEx2_results),ALDEx2_results$effect))# greater than 1 is recommended as reproducible significance cutoff
colnames(df1) <- c("seqid","ALDEx2_adj")
# ANCOMBC_result
df2 <- as.data.frame(cbind(ANCOMBC_results$res$taxon,select(ANCOMBC_results$res, starts_with("p_"), -starts_with("p_(")),select(ANCOMBC_results$res, starts_with("q_"), -starts_with("q_("))))
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
dfall <- dflist %>% reduce(full_join, by='seqid')
filename7 <- paste0(studyname,"_merged_all_result.txt")
write.table(dfall, file=filename7, quote=FALSE, sep='\t', col.names = NA)

# next we will take the combined result of that function "dfall" 
# and present results 
# split full table by Pvalue and adj Pvalue
adjustedOnly <- select(dfall, ends_with("_adj"))
PvalueOnly <- select(dfall, ends_with("_Pvalue"))
row.names(adjustedOnly) <- dfall$seqid
row.names(PvalueOnly) <- dfall$seqid
write.table(PvalueOnly, file = "PvalueOnly.txt", sep = "\t", quote = FALSE)
write.table(adjustedOnly, file = "adjustedOnly.txt", sep = "\t", quote = FALSE)
adjustedOnly <- read.table("adjustedOnly.txt", header = TRUE, row.names = 1)
PvalueOnly <- read.table("PvalueOnly.txt", header = TRUE, row.names = 1)

# heatmaps but not terribly useful in larger datasets
p2 <- filter(PvalueOnly, rowSums(is.na(PvalueOnly))<=1)
a2 <- filter(adjustedOnly, rowSums(is.na(adjustedOnly))<=1)
superheat(p2, pretty.order.rows = TRUE, pretty.order.cols = TRUE, heat.pal = c("#2D5C4F","#829D7D","#5E8CF"),heat.pal.values = c(0,0.05,1))
superheat(a2[,2:6], pretty.order.rows = TRUE, pretty.order.cols = TRUE, heat.pal = c("#2D5C4F","#829D7D","#E5E8CF"),heat.pal.values = c(0,0.05,1))
# just for any tool returning p<0.05
p2 <- filter_all(PvalueOnly, any_vars(. < 0.05))
a2 <- filter_all(adjustedOnly[,2:6], any_vars(. < 0.05))
write.table(p2, file = "PvalueOnly_Any_tool_significant_only.txt", sep = "\t", quote = FALSE)
write.table(a2, file = "adjustedOnly_Any_tool_significant_only.txt", sep = "\t", quote = FALSE)
p2 <- filter_all(PvalueOnly, any_vars(. < 0.05))
superheat(p2, pretty.order.rows = TRUE, pretty.order.cols = TRUE, heat.pal = c("#2D5C4F","#829D7D","#E5E8CF"),heat.pal.values = c(0,0.05,1))
superheat(a2[,2:6], pretty.order.rows = TRUE, pretty.order.cols = TRUE, heat.pal = c("#2D5C4F","#829D7D","#E5E8CF"),heat.pal.values = c(0,0.05,1))

# get only p value <= 0.05 more than 1 tool
PvalueOnly$count <- rowSums(PvalueOnly <= 0.05)
ConsensusPvalue <- subset(PvalueOnly, count > 1)
ConsensusPvalue <- select(ConsensusPvalue, -ends_with("count"))
write.table(ConsensusPvalue, file = "ConsensusSigPvalueOnly.txt", sep = "\t", quote = FALSE)

adjustedOnly$count <- rowSums(adjustedOnly <= 0.05)
ConsensusAdjustedOnly <- subset(adjustedOnly, count > 1)
ConsensusAdjustedOnly <- select(ConsensusAdjustedOnly, -ends_with("count"))
write.table(ConsensusAdjustedOnly, file = "ConsensusSigadjustedOnly.txt", sep = "\t", quote = FALSE)
vectorOfTables <- vector(mode = "list", length = N)
vectorOfTables[[1]] <- ConsensusPvalue
vectorOfTables[[2]] <- ConsensusAdjustedOnly
names(vectorOfTables) <- c("ConsensusPvalue","ConsensusAdjustedOnly")
return(vectorOfTables)
}

# From there we could take the significant ASVs and plot the counts per sample vs Group
both <- merge(asv2,groupings, by = "row.names")
both <- column_to_rownames(both, var = "Row.names")
melted <- melt(both)
ggplot(melted, aes(Group,log(value))) + geom_boxplot() + facet_wrap(~variable) + theme_bw() + ylab("log(count)")
