![logo](/images/Picture1.png)

**Processing of raw sequence data**

  Each Method (16S/ITS/18S) is treated differently
  1) **trim_(16S/ITS/18S).sh**
       
       runs cutadapt to remove adapters and primers
       primers are specific to your experimental design so adjust accordingly
  2) **run_(16S/ITS/18S)_dada2.R**
       runs dada2 on your samples
       
       filterAndTrim is specific to the length of your reads and type of read
       
       mergePairs overlap is specific to the variable region length
       
       taxonony references are specific to your experimental design 
       
       see: https://benjjneb.github.io/dada2/training.html for reference files
       
       recommend using references with outgroups ie ("All eukaryotes") when possible. 



**R - Analysis of dada2 output**

**You need the following files**
1) seqtab output file from dada2
2) metadata with sample names matching seqtab output **PROVIDED BY USER**
3) taxonomy file with sequences matching seqtab output and taxa output from IDtaxa in R
```
# input dada2 sequence table
t1 <- readRDS('seqtab_nochim.rds')
# input metadata
t2 <- read.table('metadata.txt', sep = '\t', comment='', head=TRUE, row.names=1, check.names = FALSE)
# input taxonomy
t3 <- readRDS("taxID.rds")

outtab <- Create_Tables(t1,t2,t3)
combined_taxa <- read.table(file = "combined_sequences_taxa.txt", sep = "\t")

taxa_out <- Make_Taxa_Tables("combined_sequences_taxa.txt")
Taxonomy_Plots(outtab$newmap)
sequence_count_table <- read.delim("sequence_process_summary.txt", row.names=1)
sequence_count_plot(sequence_count_table)
brayWmeta <- diversity(outtab$newmap,outtab$newtable)
Diversity_Plots(brayWmeta, outtab$newmap)
```
**The output from this pipeline is:**
1) combined tables for use in analysis
2) plots
 - alpha diversity
 - beta diversity
 - taxonomy plots
3) statistical analysis
4) summary plots
