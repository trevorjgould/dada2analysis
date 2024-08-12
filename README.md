![logo](/inst/logos/Picture1.png)

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

  3) **blastPlus taxonomy addition**
  		edit slurm script [blastplus_dada2.sh] and run blastplus_dada2.sh
  		
  		- This requires access to blastplus module and nt database
  		
  		then run blast_out_merge_dada2.R

**R - Analysis of dada2 output**

**You need the following files**
1) seqtab output file from dada2
2) metadata with sample names matching seqtab output **PROVIDED BY USER**
3) taxonomy file with sequences matching seqtab output and taxa output from IDtaxa in R
```
library(dada2analysis)
# load tables
inputtable <- readRDS("seqtab_nochim.rds")
taxa <- readRDS("taxa.rds")
metadata <- read.table("metadata.txt", header= TRUE, row.names = 1)

outtab <- Create_Tables(inputtable,metadata,taxa)
brayWmeta <- dada2analysis::diversity(outtab$newmap,outtab$newtable)
outlist <- rarefyTable(outtab$newtable)

# list to single objects
PCOA <- outlist$PCOA
brayWmeta <- outlist$brayWmeta
ds.mes <- outlist$ds.mes

Diversity_Plots(brayWmeta, outtab$newmap)
taxaout <- dada2analysis::Make_Taxa_Tables(outtab$combined_taxa)
Create_Taxonomy_Plot(taxout$PT,10,outtab$newmap)
Create_Taxonomy_Plot(taxout$CT,10,outtab$newmap)
Create_Taxonomy_Plot(taxout$OT,10,outtab$newmap)
Create_Taxonomy_Plot(taxout$FT,10,outtab$newmap)
Create_Taxonomy_Plot(taxout$GT,10,outtab$newmap)
Create_Taxonomy_Plot(taxout$ST,10,outtab$newmap)
```
**The output from this pipeline is:**
1) combined tables for use in analysis
2) plots
 - alpha diversity
 - beta diversity
 - taxonomy plots
3) statistical analysis
4) summary plots
