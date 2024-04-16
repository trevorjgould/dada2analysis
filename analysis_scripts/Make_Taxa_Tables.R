#' Make Taxa Tables
#' This function creates taxa tables from dada2 output
#' table of metadata processed by Make_Tables.
#'
#' @importFrom plyr ddply
#' @importFrom tibble rownames_to_column
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export

# avoiding no visible binding for global variable
domain <- phylum <- family <- genus <- species <- NULL

# reads in table from Make_Tables.R
#Make_Taxa_Tables <- function(x){
combined_taxa <- read.table("combined_sequences_taxa.txt", sep = "\t", check.names = FALSE)
#make split taxa tables
levels <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
n <- (ncol(combined_taxa) - 7)
#`%>%` <- dplyr::`%>%`
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
# remove columns of sum = 0
#KT <- KT[,colSums(KT)>0]
#PT <- PT[,colSums(PT)>0]
#CT <- CT[,colSums(CT)>0]
#OT <- OT[,colSums(OT)>0]
#FT <- FT[,colSums(FT)>0]
#GT <- GT[,colSums(GT)>0]
#ST <- ST[,colSums(ST)>0]
# remove "Unknown NA" column
KT <- KT[ , -which(names(KT) %in% c("Unknown NA"))]
PT <- PT[ , -which(names(PT) %in% c("Unknown NA"))]
CT <- CT[ , -which(names(CT) %in% c("Unknown NA"))]
OT <- OT[ , -which(names(OT) %in% c("Unknown NA"))]
FT <- FT[ , -which(names(FT) %in% c("Unknown NA"))]
GT <- GT[ , -which(names(GT) %in% c("Unknown NA"))]
ST <- ST[ , -which(names(ST) %in% c("Unknown NA"))]
KT <- KT %>% rownames_to_column("SampleID")
PT <- PT %>% rownames_to_column("SampleID")
CT <- CT %>% rownames_to_column("SampleID")
OT <- OT %>% rownames_to_column("SampleID")
FT <- FT %>% rownames_to_column("SampleID")
GT <- GT %>% rownames_to_column("SampleID")
ST <- ST %>% rownames_to_column("SampleID")
write_tsv(KT, file = "Kingdom_taxonomy.txt")
write_tsv(PT, file = "Phylum_taxonomy.txt")
write_tsv(CT, file = "Class_taxonomy.txt")
write_tsv(OT, file = "Order_taxonomy.txt")
write_tsv(FT, file = "Family_taxonomy.txt")
write_tsv(GT, file = "Genus_taxonomy.txt")
write_tsv(ST, file = "Species_taxonomy.txt")
return(list(PT=PT,CT=CT,OT=OT,FT=FT,GT=GT,ST=ST))
}
