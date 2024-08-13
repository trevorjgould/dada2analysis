#' Make Taxa Tables with Species
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
Domain_table <- combined_taxa %>% dplyr::select(Kingdom,all_of(1:n))
Phylum_table <- combined_taxa %>% dplyr::select(Phylum,all_of(1:n))
Class_table <- combined_taxa %>% dplyr::select(Class,all_of(1:n))
Order_table <- combined_taxa %>% dplyr::select(Order,all_of(1:n))
Family_table <- combined_taxa %>% dplyr::select(Family,all_of(1:n))
Genus_table <- combined_taxa %>% dplyr::select(Genus,all_of(1:n))
Species_table <- combined_taxa %>% dplyr::select(Species,all_of(1:n))
Domain <- plyr::ddply(Domain_table, "combined_taxa$Kingdom", plyr::numcolwise(sum))
Phylum <- plyr::ddply(Phylum_table, "combined_taxa$Phylum", plyr::numcolwise(sum))
Class <- plyr::ddply(Class_table, "combined_taxa$Class", plyr::numcolwise(sum))
Order <- plyr::ddply(Order_table, "combined_taxa$Order", plyr::numcolwise(sum))
Family <- plyr::ddply(Family_table, "combined_taxa$Family", plyr::numcolwise(sum))
Genus <- plyr::ddply(Genus_table, "combined_taxa$Genus", plyr::numcolwise(sum))
Species <- plyr::ddply(Species_table, "combined_taxa$Species", plyr::numcolwise(sum))
Domain = setNames(data.frame(t(Domain[,-1])), Domain[,1])
Phylum = setNames(data.frame(t(Phylum[,-1])), Phylum[,1])
Class = setNames(data.frame(t(Class[,-1])), Class[,1])
Order = setNames(data.frame(t(Order[,-1])), Order[,1])
Family = setNames(data.frame(t(Family[,-1])), Family[,1])
Genus = setNames(data.frame(t(Genus[,-1])), Genus[,1])
Species = setNames(data.frame(t(Species[,-1])), Species[,1])
# remove columns of sum = 0
#Domain <- Domain[,colSums(Domain)>0]
#Phylum <- Phylum[,colSums(Phylum)>0]
#Class <- Class[,colSums(Class)>0]
#Order <- Order[,colSums(Order)>0]
#Family <- Family[,colSums(Family)>0]
#Genus <- Genus[,colSums(Genus)>0]
#Species <- Species[,colSums(Species)>0]
# remove "Unknown NA" column
if("Unknown NA" %in% colnames(Domain)){Domain <- Domain[ , -which(names(Domain) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Phylum)){Phylum <- Phylum[ , -which(names(Phylum) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Class)){Class <- Class[ , -which(names(Class) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Order)){Order <- Order[ , -which(names(Order) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Family)){Family <- Family[ , -which(names(Family) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Genus)){Genus <- Genus[ , -which(names(Genus) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Species)){Species <- Species[ , -which(names(Species) %in% c("Unknown NA"))]}
if(NA %in% colnames(Domain)){Domain <- Domain[!is.na(names(Domain))]]}
Domain2 <- Domain %>% tibble::rownames_to_column("SampleID")
Phylum2 <- Phylum %>% tibble::rownames_to_column("SampleID")
Class2 <- Class %>% tibble::rownames_to_column("SampleID")
Order2 <- Order %>% tibble::rownames_to_column("SampleID")
Family2 <- Family %>% tibble::rownames_to_column("SampleID")
Genus2 <- Genus %>% tibble::rownames_to_column("SampleID")
Species2 <- Species %>% tibble::rownames_to_column("SampleID")
readr::write_tsv(Domain2, file = "Kingdom_taxonomy.txt")
readr::write_tsv(Phylum2, file = "Phylum_taxonomy.txt")
readr::write_tsv(Class2, file = "Class_taxonomy.txt")
readr::write_tsv(Order2, file = "Order_taxonomy.txt")
readr::write_tsv(Family2, file = "Family_taxonomy.txt")
readr::write_tsv(Genus2, file = "Genus_taxonomy.txt")
readr::write_tsv(Species2, file = "Species_taxonomy.txt")
return(list(Domain=Domain,Phylum=Phylum,Class=Class,Order=Order,Family=Family,Genus=Genus,Species=Species))
}
