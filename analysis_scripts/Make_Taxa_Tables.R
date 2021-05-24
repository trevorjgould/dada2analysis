#' Make Taxa Tables
#' This function creates taxa tables from dada2 output
#' table of metadata processed by Make_Tables.
#'
#' @importFrom plyr ddply
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export

# avoiding no visible binding for global variable
domain <- phylum <- family <- genus <- species <- NULL

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
