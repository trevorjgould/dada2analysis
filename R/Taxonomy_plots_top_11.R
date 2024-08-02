#' Taxonomy Plots
#'
#' This function creates plots that summarize prevalent taxa
#'
#' table of metadata.
#'
#' @export


# avoiding no visible binding for global variable
value <- variable <- Samples <- NULL


Taxonomy_Plots <- function(meta){
  # read in tables
  pal = c('#91bfdb','#ffffbf','#fc8d59')
  #meta <- read.table("Metadata_common.txt", sep = "\t", check.names = FALSE)
  KT <- read.table(file = "Kingdom_taxonomy.txt", sep = "\t", check.names = FALSE, header=TRUE, row.names = 1)
  PT <- read.table(file = "Phylum_taxonomy.txt", sep = "\t", check.names = FALSE, header=TRUE, row.names = 1)
  CT <- read.table(file = "Class_taxonomy.txt", sep = "\t", check.names = FALSE, header=TRUE, row.names = 1)
  OT <- read.table(file = "Order_taxonomy.txt", sep = "\t", check.names = FALSE, header=TRUE, row.names = 1)
  FT <- read.table(file = "Family_taxonomy.txt", sep = "\t", check.names = FALSE, header=TRUE, row.names = 1)
  GT <- read.table(file = "Genus_taxonomy.txt", sep = "\t", check.names = FALSE, header=TRUE, row.names = 1)
  ST <- read.table(file = "Species_taxonomy.txt", sep = "\t", check.names = FALSE, header=TRUE, row.names = 1)
  #meta <- outtab$newmap[-which(outtab$newmap$Eligibility_group == ""), ]
  #meta <- outtab$newmap[-which(outtab$newmap$V1_ACV2SPIKE_Result == ""), ]

  # Here we are taking the taxonomy tables created above
  # taxa with less than a sum of 0.1 total proportion over all samples are merged
  # into a column of OTHER with the NA column if it is present.
  # NEXT IS TO DYNAMICALLY SET THE 0.1 TO LEAVE 20 TAXA TOTAL WITH OTHER AS REST.
  taxa_list <- list(KT,PT,CT,OT,FT,GT,ST)
  taxa_names <- list("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  for (x in 1:7){
    XT_other = 0
    xt=0
    XTtax=0
    name_label <- taxa_names[x]
    xt=as.data.frame(taxa_list[x])
    #rownames(xt) <- xt[,1]
    #xt <- xt[,-c(1)]
    #xt=xt[(rowSums(xt)>0),]
    
    # fix rownames . to -
    row.names(meta) <- gsub("-", "\\.", row.names(meta))
    # specific to leading X in rownames
    # row.names(xt) <- gsub("X27", "27", row.names(xt))
    #find common names
    common <- intersect(rownames(meta),rownames(xt))
    # get just the overlapping samples
    meta <- meta[common,, drop = FALSE]
    xt <- xt[common,, drop = FALSE]

    XTtax = sweep(xt, 1, rowSums(xt),'/')
    # KEEP TOP TAXA AND SORT REST INTO OTHER CATEGORY AND PLOT
    # get top 15 taxa
    keepers <- row.names(as.data.frame(head(sort(colSums(XTtax), decreasing = TRUE), n=15)))
    XT_other <- XTtax[,keepers]
    # sum taxa columns smaller then top 12
    if (ncol(XTtax) == 16){XT_other = XTtax
    } else if (ncol(XTtax) > 15){
      XT_other$other <- rowSums(XTtax[,-which(names(XTtax) %in% keepers)])
      #
      if("NA." %in% colnames(XT_other))
      {
        cat("Merging NA column with Other!\n");
        #16S
        XT_other$Other <- XT_other$'NA.' + XT_other$other
        # FUNGI
        #XT_other$Other <- XT_other$'Unknown Fungi' + XT_other$other
        # remove columns NA and other
        #16S
        XT_other <- XT_other[ , -which(names(XT_other) %in% c("NA.","other"))]
        # FUNGI
        #XT_other <- XT_other[ , -which(names(XT_other) %in% c("Unknown Fungi","other"))]
      }
      b <- sum(nrow(XT_other))
      a <- as.integer(sum(rowSums(XT_other)))
      
      if(a==b){
        cat("taxa are looking good\n")
      }
      if(a!=b){
        cat("taxa do not sum to 100. there is something wrong\n")
      }
    } else {XT_other = XTtax}
  
    
    both <- cbind(XT_other,meta)
    both$Samples <- row.names(meta)
    melted <- reshape2::melt(both, id.vars = c(colnames(meta),"Samples"))
    melted$variable <- gsub("[a-z]__", "", melted$variable)
    melted$variable <- gsub("\\.", " ", melted$variable)
    filename <- paste0(name_label,"_taxonomy_other.png")
    filename2 <- paste0(name_label,"_taxonomy_other.txt")
    #p <- ggplot2::ggplot(melted, ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + facet_grid(~as.factor(V1), scales = "free") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + ggplot2::scale_fill_discrete(name = name_label)
    #ggplot2::ggsave(p, file = filename, dpi  = 800, width = 10, height = 8, units = "in")
    write.table(both, file = filename2, sep = "\t", quote = FALSE)
    #set colors
    library(RColorBrewer)
    colourCount = 16
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    #Sites
    melted$variable = with(melted, reorder(variable,value,mean))
    melted$SampleID <- rownames(melted)
    p2 <- ggplot(data = melted[!is.na(melted$SampleID),], ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + scale_x_discrete(expand = c(0, 0.5)) + facet_grid(~Treatment, scales = "free") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) + scale_fill_manual(values = getPalette(colourCount),name = name_label) + ggtitle("SampleID")
    filename3 <- paste0(name_label,"Sample_taxonomy_Treatment.png")
    ggplot2::ggsave(p2, file = filename3, dpi  = 800, width = 8, height = 6, units = "in")
    
    #p3 <- ggplot(data = melted[!is.na(melted$SampleID),], ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + scale_x_discrete(expand = c(0, 0.5)) + facet_grid(~V1, scales = "free") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) + scale_fill_manual(values = getPalette(colourCount),name = name_label)
    #filename4 <- paste0(name_label,"Sample_taxonomy_V1.png")
    #ggplot2::ggsave(p3, file = filename4, dpi  = 800, width = 8, height = 6, units = "in")
   
    #p4 <- ggplot(data = melted[!is.na(melted$SampleID),], ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + scale_x_discrete(expand = c(0, 0.5)) + facet_grid(~V2, scales = "free") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) + scale_fill_manual(values = getPalette(colourCount),name = name_label)
    #filename5 <- paste0(name_label,"Sample_taxonomy_V2.png")
    #ggplot2::ggsave(p4, file = filename5, dpi  = 800, width = 8, height = 6, units = "in")
    
    #p4 <- ggplot(data = melted[!is.na(melted$SampleID),], ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + scale_x_discrete(expand = c(0, 0.5)) + facet_grid(~V3, scales = "free") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) + scale_fill_manual(values = getPalette(colourCount),name = name_label)
    #filename5 <- paste0(name_label,"Sample_taxonomy_V3.png")
    #ggplot2::ggsave(p4, file = filename5, dpi  = 800, width = 8, height = 6, units = "in")
    # heatmap
    #melted$variable = with(melted, reorder(variable,value,mean))
    #p2 <- ggplot(melted, aes(Samples,variable, fill=value)) + geom_tile(color="white", size=0.1) + facet_grid(~Sites, scales = "free_y") + xlab("Day") + ylab(name_label) + theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colours = pal,name = "Proportion") + ggtitle("Sites")
    
  }
}

