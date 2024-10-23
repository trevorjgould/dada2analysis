library(dada2analysis)
# load tables
inputtable <- readRDS("seqtab_nochim.rds")
taxa <- readRDS("taxa.rds")
metadata <- read.table("metadata.txt", header= TRUE, row.names = 1)

outtab <- Create_Tables(inputtable,metadata,taxa)
brayWmeta <- diversity(outtab$newmap,outtab$newtable)
outlist <- rarefyTable(outtab$newtable)

# list to single objects
PCOA <- outlist$PCOA
brayWmeta <- outlist$brayWmeta
ds.mes <- outlist$ds.mes
combined_taxa = outtab$combined_taxa
Diversity_Plots(brayWmeta, outtab$newmap)
taxaout <- Make_Taxa_Tables(combined_taxa)
Create_Taxonomy_Plot(taxaout$Phylum,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Class,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Order,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Family,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Genus,10,outtab$newmap)
Create_Taxonomy_Plot(taxaout$Species,10,outtab$newmap)
