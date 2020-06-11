This is the analysis pipeline for use with dada2 processed output 
you need the following. 
1) seqtab_output from dada2
2) metadata with sample names matching seqtab output
3) taxonomy file with sequences matching seqtab output and taxa output from IDtaxa in R

The output from this pipeline is
1) combined tables for use in analysis
2) plots
 - alpha diversity
 - beta diversity
 - taxonomy plots
3) statistical analysis
4) summary plots

Processing of 16S sequences
1) dada2_16S_processing.pbs
calls:
2) dada2_version2.R

Run Dada2
input
- sequence data
output
- seqtab_nochim.rds
- taxID.rds
User Provided
- metadata.txt

1) Make_Tables.Rmd
input
- seqtab_nochim.rds
- metadata.txt
- taxID.rds
output
- Sequence_table_common.rds
- Metadata_common.txt
- combined_taxa.txt

2) Make_Taxa_Tables.Rmd
input: 	
- Metadata_common.txt
- combined_taxa.txt
output:	
- Kingdom_taxonomy.txt
- Phylum_taxonomy.txt
- Class_taxonomy.txt
- Order_taxonomy.txt
- Family_taxonomy.txt
- Genus_taxonomy.txt

3) Taxonomy_Plots.Rmd
input:
- Kingdom_taxonomy.txt
- Phylum_taxonomy.txt
- Class_taxonomy.txt
- Order_taxonomy.txt
- Family_taxonomy.txt
- Genus_taxonomy.txt
- proportional_diversity_stats.txt
output:	
- [Level]_taxonomy_other.png