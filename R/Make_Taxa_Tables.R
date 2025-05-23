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
Make_Taxa_Tables <- function(combined_taxa){
Kingdom <- Phylum <- Class <- Order <- Family <- Genus <- NULL
#combined_taxa <- read.table("combined_sequences_taxa.txt", sep = "\t", check.names = FALSE)
#make split taxa tables
levels <- c("Kingdom","Phylum","Class","Order","Family","Genus")
n <- (ncol(combined_taxa) - 6)
combined_taxa[,1:n] <- sapply(combined_taxa[,1:n], as.numeric)
#
Domain_table <- combined_taxa %>% dplyr::select(Kingdom,all_of(1:n))
Phylum_table <- combined_taxa %>% dplyr::select(Phylum,all_of(1:n))
Class_table <- combined_taxa %>% dplyr::select(Class,all_of(1:n))
Order_table <- combined_taxa %>% dplyr::select(Order,all_of(1:n))
Family_table <- combined_taxa %>% dplyr::select(Family,all_of(1:n))
Genus_table <- combined_taxa %>% dplyr::select(Genus,all_of(1:n))
Domain <- plyr::ddply(Domain_table, "combined_taxa$Kingdom", plyr::numcolwise(sum))
Phylum <- plyr::ddply(Phylum_table, "combined_taxa$Phylum", plyr::numcolwise(sum))
Class <- plyr::ddply(Class_table, "combined_taxa$Class", plyr::numcolwise(sum))
Order <- plyr::ddply(Order_table, "combined_taxa$Order", plyr::numcolwise(sum))
Family <- plyr::ddply(Family_table, "combined_taxa$Family", plyr::numcolwise(sum))
Genus <- plyr::ddply(Genus_table, "combined_taxa$Genus", plyr::numcolwise(sum))
Domain = setNames(data.frame(t(Domain[,-1])), Domain[,1])
Phylum = setNames(data.frame(t(Phylum[,-1])), Phylum[,1])
Class = setNames(data.frame(t(Class[,-1])), Class[,1])
Order = setNames(data.frame(t(Order[,-1])), Order[,1])
Family = setNames(data.frame(t(Family[,-1])), Family[,1])
Genus = setNames(data.frame(t(Genus[,-1])), Genus[,1])
# remove columns of sum = 0
#Domain <- Domain[,colSums(Domain)>0]
#Phylum <- Phylum[,colSums(Phylum)>0]
#Class <- Class[,colSums(Class)>0]
#Order <- Order[,colSums(Order)>0]
#Family <- Family[,colSums(Family)>0]
#Genus <- Genus[,colSums(Genus)>0]
# remove "Unknown NA" column
if("Unknown NA" %in% colnames(Domain)){Domain <- Domain[ , -which(names(Domain) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Phylum)){Phylum <- Phylum[ , -which(names(Phylum) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Class)){Class <- Class[ , -which(names(Class) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Order)){Order <- Order[ , -which(names(Order) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Family)){Family <- Family[ , -which(names(Family) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Genus)){Genus <- Genus[ , -which(names(Genus) %in% c("Unknown NA"))]}
if(NA %in% colnames(Domain)){Domain <- Domain[!is.na(names(Domain))]}
Domain = Domain[,order(colSums(-Domain))]
Phylum = Phylum[,order(colSums(-Phylum))]
Class = Class[,order(colSums(-Class))]
Order = Order[,order(colSums(-Order))]
Family = Family[,order(colSums(-Family))]
Genus = Genus[,order(colSums(-Genus))]
Domain2 <- Domain %>% tibble::rownames_to_column("SampleID")
Phylum2 <- Phylum %>% tibble::rownames_to_column("SampleID")
Class2 <- Class %>% tibble::rownames_to_column("SampleID")
Order2 <- Order %>% tibble::rownames_to_column("SampleID")
Family2 <- Family %>% tibble::rownames_to_column("SampleID")
Genus2 <- Genus %>% tibble::rownames_to_column("SampleID")
readr::write_tsv(Domain2, file = "Kingdom_taxonomy.txt")
readr::write_tsv(Phylum2, file = "Phylum_taxonomy.txt")
readr::write_tsv(Class2, file = "Class_taxonomy.txt")
readr::write_tsv(Order2, file = "Order_taxonomy.txt")
readr::write_tsv(Family2, file = "Family_taxonomy.txt")
readr::write_tsv(Genus2, file = "Genus_taxonomy.txt")
return(list(Domain=Domain,Phylum=Phylum,Class=Class,Order=Order,Family=Family,Genus=Genus))
}

Make_Taxa_Tables_18S <- function(combined_taxa){
Domain <- Supergroup <- Division <- Subdivision <- Class <- Order <- Family <- Genus <- Species <- NULL
#combined_taxa <- read.table("combined_sequences_taxa.txt", sep = "\t", check.names = FALSE)
#make split taxa tables
levels <- c("Domain","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species")
n <- (ncol(combined_taxa) - 9)
combined_taxa[,1:n] <- sapply(combined_taxa[,1:n], as.numeric)
#
Domain_table <- combined_taxa %>% dplyr::select(Domain,all_of(1:n))
Supergroup_table <- combined_taxa %>% dplyr::select(Supergroup,all_of(1:n))
Division_table <- combined_taxa %>% dplyr::select(Division,all_of(1:n))
Subdivision_table <- combined_taxa %>% dplyr::select(Subdivision,all_of(1:n))
Class_table <- combined_taxa %>% dplyr::select(Class,all_of(1:n))
Order_table <- combined_taxa %>% dplyr::select(Order,all_of(1:n))
Family_table <- combined_taxa %>% dplyr::select(Family,all_of(1:n))
Genus_table <- combined_taxa %>% dplyr::select(Genus,all_of(1:n))
Domain <- plyr::ddply(Domain_table, "combined_taxa$Domain", plyr::numcolwise(sum))
Supergroup <- plyr::ddply(Supergroup_table, "combined_taxa$Supergroup", plyr::numcolwise(sum))
Division <- plyr::ddply(Division_table, "combined_taxa$Division", plyr::numcolwise(sum))
Subdivision <- plyr::ddply(Subdivision_table, "combined_taxa$Subdivision", plyr::numcolwise(sum))
Class <- plyr::ddply(Class_table, "combined_taxa$Class", plyr::numcolwise(sum))
Order <- plyr::ddply(Order_table, "combined_taxa$Order", plyr::numcolwise(sum))
Family <- plyr::ddply(Family_table, "combined_taxa$Family", plyr::numcolwise(sum))
Genus <- plyr::ddply(Genus_table, "combined_taxa$Genus", plyr::numcolwise(sum))
Domain = setNames(data.frame(t(Domain[,-1])), Domain[,1])
Supergroup = setNames(data.frame(t(Supergroup[,-1])), Supergroup[,1])
Division = setNames(data.frame(t(Division[,-1])), Division[,1])
Subdivision = setNames(data.frame(t(Subdivision[,-1])), Subdivision[,1])
Class = setNames(data.frame(t(Class[,-1])), Class[,1])
Order = setNames(data.frame(t(Order[,-1])), Order[,1])
Family = setNames(data.frame(t(Family[,-1])), Family[,1])
Genus = setNames(data.frame(t(Genus[,-1])), Genus[,1])
# remove columns of sum = 0
#Domain <- Domain[,colSums(Domain)>0]
#Supergroup <- Supergroup[,colSums(Supergroup)>0]
#Class <- Class[,colSums(Class)>0]
#Order <- Order[,colSums(Order)>0]
#Family <- Family[,colSums(Family)>0]
#Genus <- Genus[,colSums(Genus)>0]
# remove "Unknown NA" column
if("Unknown NA" %in% colnames(Domain)){Domain <- Domain[ , -which(names(Domain) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Supergroup)){Supergroup <- Supergroup[ , -which(names(Supergroup) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Division)){Division <- Division[ , -which(names(Division) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Subdivision)){Subdivision <- Subdivision[ , -which(names(Subdivision) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Class)){Class <- Class[ , -which(names(Class) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Order)){Order <- Order[ , -which(names(Order) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Family)){Family <- Family[ , -which(names(Family) %in% c("Unknown NA"))]}
if("Unknown NA" %in% colnames(Genus)){Genus <- Genus[ , -which(names(Genus) %in% c("Unknown NA"))]}
if(NA %in% colnames(Domain)){Domain <- Domain[!is.na(names(Domain))]}
Domain = Domain[,order(colSums(-Domain))]
Supergroup = Supergroup[,order(colSums(-Supergroup))]
Division = Division[,order(colSums(-Division))]
Subdivision = Subdivision[,order(colSums(-Subdivision))]
Class = Class[,order(colSums(-Class))]
Order = Order[,order(colSums(-Order))]
Family = Family[,order(colSums(-Family))]
Genus = Genus[,order(colSums(-Genus))]
Domain2 <- Domain %>% tibble::rownames_to_column("SampleID")
Supergroup2 <- Supergroup %>% tibble::rownames_to_column("SampleID")
Division2 <- Division %>% tibble::rownames_to_column("SampleID")
Subdivision2 <- Subdivision %>% tibble::rownames_to_column("SampleID")
Class2 <- Class %>% tibble::rownames_to_column("SampleID")
Order2 <- Order %>% tibble::rownames_to_column("SampleID")
Family2 <- Family %>% tibble::rownames_to_column("SampleID")
Genus2 <- Genus %>% tibble::rownames_to_column("SampleID")
readr::write_tsv(Domain2, file = "Domain_taxonomy.txt")
readr::write_tsv(Supergroup2, file = "Supergroup_taxonomy.txt")
readr::write_tsv(Division2, file = "Division_taxonomy.txt")
readr::write_tsv(Subdivision2, file = "Subdivision_taxonomy.txt")
readr::write_tsv(Class2, file = "Class_taxonomy.txt")
readr::write_tsv(Order2, file = "Order_taxonomy.txt")
readr::write_tsv(Family2, file = "Family_taxonomy.txt")
readr::write_tsv(Genus2, file = "Genus_taxonomy.txt")
return(list(Domain=Domain,Supergroup=Supergroup,Division=Division,Subdivision=Subdivision,Class=Class,Order=Order,Family=Family,Genus=Genus))
}
