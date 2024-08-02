#' Taxonomy Plots
#'
#' This function creates plots that summarize prevalent taxa
#'
#' table of metadata.
#'
#' @export


# avoiding no visible binding for global variable
value <- variable <- Samples <- NULL


Taxonomy_Plots <- function(h){
# read in tables
pal = c('#91bfdb','#ffffbf','#fc8d59')
meta <- outtab$newmap
#KT <- read.table(file = "Kingdom_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
PT <- read.table(file = "Phylum_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
CT <- read.table(file = "Class_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
OT <- read.table(file = "Order_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
FT <- read.table(file = "Family_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
GT <- read.table(file = "Genus_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
ST <- read.table(file = "Species_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
#Group = meta %>% select(h)
Gname <- colnames(Group)

# Here we are taking the taxonomy tables created above
# taxa with less than a sum of 0.1 total proportion over all samples are merged
# into a column of OTHER with the NA column if it is present.
# NEXT IS TO DYNAMICALLY SET THE 0.1 TO LEAVE 20 TAXA TOTAL WITH OTHER AS REST.
taxa_list <- list(KT,PT,CT,OT,FT,GT,ST)
taxa_names <- list("Kingdom","Phylum","Class","Order","Family","Genus","Species")
for (x in 3:7){
XT_other = 0
xt=0
XTtax=0
name_label <- taxa_names[x]
xt=as.data.frame(taxa_list[x])
xt=xt[(rowSums(xt)>0),]

# fix rownames . to -
row.names(meta) <- gsub("-", "\\.", row.names(meta))

#find common names
common <- intersect(rownames(meta),rownames(xt))
# get just the overlapping samples
meta <- meta[common,, drop = FALSE]
xt <- xt[common,, drop = FALSE]


XTtax = sweep(xt, 1, rowSums(xt),'/')
# KEEP TOP TAXA AND SORT REST INTO OTHER CATEGORY AND PLOT
# get taxa that are present > 0.1
XT_other <- XTtax[,(colSums(XTtax)>=0.1)]
# sum taxa less then 0.1

if (table(colSums(XTtax)<0.1)["TRUE"] > 1){
XT_other$other <- rowSums(XTtax[,(colSums(XTtax)<0.1)])

#
if("NA." %in% colnames(XT_other))
{
cat("Merging NA column with Other!\n");
XT_other$Other <- XT_other$'NA.' + XT_other$other
# remove columns NA and other
XT_other <- XT_other[ , -which(names(XT_other) %in% c("NA.","other"))]
}
b <- sum(nrow(XT_other))
a <- as.integer(sum(rowSums(XT_other)))

if(a==b){
cat("taxa are looking good\n")
}
if(a!=b){
 cat("taxa do not sum to 100. there is something wrong\n")
}
}


both <- cbind(XT_other,meta)
both$Samples <- row.names(meta)
melted <- reshape2::melt(both, id.vars = c(colnames(meta),"Samples"))
melted$variable <- gsub("[a-z]__", "", melted$variable)
melted$variable <- gsub("\\.", " ", melted$variable)
filename <- paste0(name_label,"_taxonomy_",Gname,".png")
filename2 <- paste0(name_label,"_taxonomy_table_",Gname,".txt")
p <- ggplot2::ggplot(melted, ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + facet_grid(~Groups, scales = "free") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x=element_blank(),axis.ticks=element_blank()) + ggplot2::scale_fill_discrete(name = name_label)
ggplot2::ggsave(p, file = filename, dpi  = 800, width = 14, height = 8, units = "in")
write.table(both, file = filename2, sep = "\t", quote = FALSE)

melted$variable = with(melted, reorder(variable,value,mean))
p2 <- ggplot(melted, aes(Samples,variable, fill=value)) + geom_tile(color="white", size=0.1) + facet_grid(~Groups, scales = "free") + ylab(name_label) + theme(axis.text.y = element_text(size = 8),axis.text.x=element_blank(),axis.ticks=element_blank()) + scale_fill_gradientn(colours = pal,name = "Proportion")
filename3 <- paste0(name_label,"_taxonomy_heat.png")
ggplot2::ggsave(p2, file = filename3, dpi  = 800, width = 6, height = 10, units = "in")
}
}

Taxonomy_Plots("Group")

Taxonomy_Plots_Top10 <- function(t,meta){
  # read in tables
  pal = c('#91bfdb','#ffffbf','#fc8d59')
  #meta <- read.table("Metadata_common.txt", sep = "\t", check.names = FALSE)
  #KT <- read.table(file = "Kingdom_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
  PT <- read.table(file = "Phylum_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
  CT <- read.table(file = "Class_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
  OT <- read.table(file = "Order_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
  FT <- read.table(file = "Family_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
  GT <- read.table(file = "Genus_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
  ST <- read.table(file = "Species_taxonomy.txt", sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)

  # Here we are taking the taxonomy tables created above
  # taxa with less than a sum of 0.1 total proportion over all samples are merged
  # into a column of OTHER with the NA column if it is present.
  # NEXT IS TO DYNAMICALLY SET THE 0.1 TO LEAVE 20 TAXA TOTAL WITH OTHER AS REST.
  taxa_list <- list(KT,PT,CT,OT,FT,GT,ST)
  taxa_names <- list("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  for (x in 3:7){
    XT_other = 0
    xt=0
    XTtax=0
    name_label <- taxa_names[x]
    xt=as.data.frame(taxa_list[x])
    xt=xt[(rowSums(xt)>0),]

    # fix rownames . to -
    row.names(meta) <- gsub("-", "\\.", row.names(meta))

    #find common names
    common <- intersect(rownames(meta),rownames(xt))
    # get just the overlapping samples
    meta <- meta[common,, drop = FALSE]
    xt <- xt[common,, drop = FALSE]


    XTtax = sweep(xt, 1, rowSums(xt),'/')
    # KEEP TOP TAXA AND SORT REST INTO OTHER CATEGORY AND PLOT
    # get taxa that are present > 0.1
    if (ncol(XTtax) > t){
      topt <- rownames(as.data.frame(head(sort(colSums(XTtax), decreasing = TRUE), n=t)))
      nottopt <- rownames(as.data.frame(sort(colSums(XTtax), decreasing = TRUE)))
      nottopt <- nottopt[-c(1:t)]
      XT_other <- XTtax[,topt]
      XT_other$other <- rowSums(XTtax[,nottopt])

      #
      if("NA." %in% colnames(XT_other))
      {
        cat("Merging NA column with Other!\n");
        XT_other$Other <- XT_other$'NA.' + XT_other$other
        # remove columns NA and other
        XT_other <- XT_other[ , -which(names(XT_other) %in% c("NA.","other"))]
      }
      b <- sum(nrow(XT_other))
      a <- as.integer(sum(rowSums(XT_other)))

      if(a==b){
        cat("taxa are looking good\n")
      }
      if(a!=b){
        cat("taxa do not sum to 100. there is something wrong\n")
      }
    }


    both <- cbind(XT_other,meta)
    both$Samples <- row.names(meta)
    melted <- reshape2::melt(both, id.vars = c(colnames(meta),"Samples"))
    melted$variable <- gsub("[a-z]__", "", melted$variable)
    melted$variable <- gsub("\\.", " ", melted$variable)
    filename <- paste0(name_label,"_taxonomy_top",t,"_Type_boxplot.png")
    filename2 <- paste0(name_label,"_taxonomy_top",t,".txt")
    titlename <- paste0("Top_",t,"_",name_label)
    forcolors <- c(rep(c("#899DA4", "#C93312"), times = 6))
    color4 = c("#BB3717","#54969A","#B47F6D","#0F5372")
    color3 = c("#A9502D","#678076","#C7723A")
    p <- ggplot2::ggplot(melted, ggplot2::aes(Group, value, fill = Group)) + ggplot2::geom_boxplot()+ ggplot2::ylab("Proportion") + ggplot2::scale_fill_manual(values = color4) + facet_grid(~variable, scales = "free", switch="both")+ theme_bw()+theme(strip.text.y.left = element_text(angle = 0),axis.text.x=element_blank(),axis.ticks=element_blank(),strip.text=element_text(margin=margin()),panel.spacing=unit(0, "lines"),legend.position="bottom") + ggtitle(titlename)
    ggplot2::ggsave(p, file = filename, dpi  = 800, width = 14, height = 8, units = "in")
    write.table(both, file = filename2, sep = "\t", quote = FALSE)
       }
}
Taxonomy_Plots_Top10(10,outtab$newmap)
