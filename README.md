![logo](/inst/logos/Picture1.png)
**Install**
devtools::install_github("trevorjgould/dada2analysis")

**You need the following files**
1) seqtab output file from dada2
2) metadata with sample names matching seqtab output **PROVIDED BY USER**
3) taxonomy file with sequences matching seqtab output and taxa output from IDtaxa in R
```
library(dada2analysis)
# LOAD TABLES
# samples as rows, ASVs as headers
inputtable <- readRDS("seqtab_nochim.rds")
# ASV sequences as rows, taxonomy levels as headers
taxa <- readRDS("taxa.rds")
# samples as rows, metadata variables as headers
metadata <- read.table("metadata.txt", header= TRUE, row.names = 1)

outtab <- Create_Tables(inputtable,metadata,taxa)
brayWmeta <- diversity(outtab$newmap,outtab$newtable)
outlist <- rarefyTable(outtab$newtable)

# list object to single objects
PCOA <- outlist$PCOA
brayWmeta <- outlist$brayWmeta
ds.mes <- outlist$ds.mes
combined_taxa <-  outtab$combined_taxa

# diversity plots
Diversity_Plots(brayWmeta, outtab$newmap)

# taxonomy tables
taxaout <- Make_Taxa_Tables(combined_taxa)

#taxonomy plots
Create_Taxonomy_Plot(taxaout$Phylum,2,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Class,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Order,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Family,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Genus,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Species,10,outtab$newmap)

Create_Taxonomy_Plot_Facet(taxaout$Phylum,10,outtab$newmap)
Create_Taxonomy_Plot_Facet(taxaout$Class,10,outtab$newmap)
Create_Taxonomy_Plot_Facet(taxaout$Order,10,outtab$newmap)
Create_Taxonomy_Plot_Facet(taxaout$Family,10,outtab$newmap)
Create_Taxonomy_Plot_Facet(taxaout$Genus,10,outtab$newmap)
Create_Taxonomy_Plot_Facet(taxaout$Species,10,outtab$newmap)
```
**The output from this pipeline is:**
1) combined tables for use in analysis
2) plots
 - alpha diversity
 - beta diversity
 - taxonomy plots
3) statistical analysis
4) summary plots
