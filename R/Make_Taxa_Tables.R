#' Make Taxa Tables
#' This function creates taxa tables from dada2 output
#' table of metadata processed by Make_Tables.
#'
#' @param combined_taxa combined table of seqtab and taxa
#' @importFrom plyr ddply
#' @importFrom tibble rownames_to_column
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export

# reads in table from Make_Tables.R
Make_Taxa_Tables_with_Species <- function(combined_taxa){
Kingdom <- Phylum <- Class <- Order <- Family <- Genus <- Species <- NULL
#combined_taxa <- read.table("combined_sequences_taxa.txt", sep = "\t", check.names = FALSE)
#make split taxa tables
levels <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
n <- (ncol(combined_taxa) - 7)
combined_taxa[,1:n] <- sapply(combined_taxa[,1:n], as.numeric)
#
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
Domain <- KT[ , -which(names(KT) %in% c("Unknown NA"))]
Phylum <- PT[ , -which(names(PT) %in% c("Unknown NA"))]
Class <- CT[ , -which(names(CT) %in% c("Unknown NA"))]
Order <- OT[ , -which(names(OT) %in% c("Unknown NA"))]
Family <- FT[ , -which(names(FT) %in% c("Unknown NA"))]
Genus <- GT[ , -which(names(GT) %in% c("Unknown NA"))]
Species <- ST[ , -which(names(ST) %in% c("Unknown NA"))]
KT2 <- KT %>% tibble::rownames_to_column("SampleID")
PT2 <- PT %>% tibble::rownames_to_column("SampleID")
CT2 <- CT %>% tibble::rownames_to_column("SampleID")
OT2 <- OT %>% tibble::rownames_to_column("SampleID")
FT2 <- FT %>% tibble::rownames_to_column("SampleID")
GT2 <- GT %>% tibble::rownames_to_column("SampleID")
ST2 <- ST %>% tibble::rownames_to_column("SampleID")
readr::write_tsv(KT2, file = "Kingdom_taxonomy.txt")
readr::write_tsv(PT2, file = "Phylum_taxonomy.txt")
readr::write_tsv(CT2, file = "Class_taxonomy.txt")
readr::write_tsv(OT2, file = "Order_taxonomy.txt")
readr::write_tsv(FT2, file = "Family_taxonomy.txt")
readr::write_tsv(GT2, file = "Genus_taxonomy.txt")
readr::write_tsv(ST2, file = "Species_taxonomy.txt")
return(list(PT=PT,CT=CT,OT=OT,FT=FT,GT=GT,ST=ST))
}

# reads in table from Make_Tables.R
Make_Taxa_Tables <- function(combined_taxa){
Kingdom <- Phylum <- Class <- Order <- Family <- Genus <- NULL
#combined_taxa <- read.table("combined_sequences_taxa.txt", sep = "\t", check.names = FALSE)
#make split taxa tables
levels <- c("Kingdom","Phylum","Class","Order","Family","Genus")
n <- (ncol(combined_taxa) - 6)
combined_taxa[,1:n] <- sapply(combined_taxa[,1:n], as.numeric)
#
Domain_table <- combined_taxa %>% dplyr::select(Kingdom,1:n)
Phylum_table <- combined_taxa %>% dplyr::select(Phylum,1:n)
Class_table <- combined_taxa %>% dplyr::select(Class,1:n)
Order_table <- combined_taxa %>% dplyr::select(Order,1:n)
Family_table <- combined_taxa %>% dplyr::select(Family,1:n)
Genus_table <- combined_taxa %>% dplyr::select(Genus,1:n)
KT <- plyr::ddply(Domain_table, "combined_taxa$Kingdom", plyr::numcolwise(sum))
PT <- plyr::ddply(Phylum_table, "combined_taxa$Phylum", plyr::numcolwise(sum))
CT <- plyr::ddply(Class_table, "combined_taxa$Class", plyr::numcolwise(sum))
OT <- plyr::ddply(Order_table, "combined_taxa$Order", plyr::numcolwise(sum))
FT <- plyr::ddply(Family_table, "combined_taxa$Family", plyr::numcolwise(sum))
GT <- plyr::ddply(Genus_table, "combined_taxa$Genus", plyr::numcolwise(sum))
KT = setNames(data.frame(t(KT[,-1])), KT[,1])
PT = setNames(data.frame(t(PT[,-1])), PT[,1])
CT = setNames(data.frame(t(CT[,-1])), CT[,1])
OT = setNames(data.frame(t(OT[,-1])), OT[,1])
FT = setNames(data.frame(t(FT[,-1])), FT[,1])
GT = setNames(data.frame(t(GT[,-1])), GT[,1])
# remove columns of sum = 0
#KT <- KT[,colSums(KT)>0]
#PT <- PT[,colSums(PT)>0]
#CT <- CT[,colSums(CT)>0]
#OT <- OT[,colSums(OT)>0]
#FT <- FT[,colSums(FT)>0]
#GT <- GT[,colSums(GT)>0]
# remove "Unknown NA" column
Domain <- KT[ , -which(names(KT) %in% c("Unknown NA"))]
Phylum <- PT[ , -which(names(PT) %in% c("Unknown NA"))]
Class <- CT[ , -which(names(CT) %in% c("Unknown NA"))]
Order <- OT[ , -which(names(OT) %in% c("Unknown NA"))]
Family <- FT[ , -which(names(FT) %in% c("Unknown NA"))]
Genus <- GT[ , -which(names(GT) %in% c("Unknown NA"))]
KT2 <- KT %>% tibble::rownames_to_column("SampleID")
PT2 <- PT %>% tibble::rownames_to_column("SampleID")
CT2 <- CT %>% tibble::rownames_to_column("SampleID")
OT2 <- OT %>% tibble::rownames_to_column("SampleID")
FT2 <- FT %>% tibble::rownames_to_column("SampleID")
GT2 <- GT %>% tibble::rownames_to_column("SampleID")
readr::write_tsv(KT2, file = "Kingdom_taxonomy.txt")
readr::write_tsv(PT2, file = "Phylum_taxonomy.txt")
readr::write_tsv(CT2, file = "Class_taxonomy.txt")
readr::write_tsv(OT2, file = "Order_taxonomy.txt")
readr::write_tsv(FT2, file = "Family_taxonomy.txt")
readr::write_tsv(GT2, file = "Genus_taxonomy.txt")
return(list(PT=PT,CT=CT,OT=OT,FT=FT,GT=GT))
}
