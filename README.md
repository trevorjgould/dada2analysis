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