**The output from this pipeline is:**
1) combined tables for use in analysis
2) plots
 - alpha diversity
 - beta diversity
 - taxonomy plots
3) statistical analysis
4) summary plots

**Processing of 16S sequences**
  - dada2_16S_processing.pbs > dada2_version2.R

- User Provided
  - metadata.txt

**Run Dada2**
input | output
--------- | ---------
sequence data | seqtab_nochim.rds
 . | taxID.rds


**You need the following files**
1) seqtab output file from dada2
2) metadata with sample names matching seqtab output **PROVIDED BY USER**
3) taxonomy file with sequences matching seqtab output and taxa output from IDtaxa in R

**Make_Tables.Rmd**
input | output
--------- | ---------
seqtab_nochim.rds | Sequence_table_common.rds
metadata.txt | Metadata_common.txt
taxID.rds | combined_taxa.txt


**Make_Taxa_Tables.Rmd**
input | output
--------- | ---------
Metadata_common.txt | Kingdom_taxonomy.txt
combined_taxa.txt | Phylum_taxonomy.txt
.  | Class_taxonomy.txt
.  | Order_taxonomy.txt
.  | Family_taxonomy.txt
.  | Genus_taxonomy.txt

**Taxonomy_Plots.Rmd**
input | output
--------- | ---------
Kingdom_taxonomy.txt | Kingdom_taxonomy_other.png
Phylum_taxonomy.txt | Phylum_taxonomy_other.png
Class_taxonomy.txt | Class_taxonomy_other.png
Order_taxonomy.txt | Order_taxonomy_other.png
Family_taxonomy.txt | Family_taxonomy_other.png
Genus_taxonomy.txt | Genus_taxonomy_other.png
Metadata_common.txt | 