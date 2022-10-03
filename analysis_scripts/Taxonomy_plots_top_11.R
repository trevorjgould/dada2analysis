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
  # KT <- read.table(file = "Kingdom_taxonomy.txt", sep = "\t", check.names = FALSE)
  PT <- read.table(file = "Phylum_taxonomy.txt", sep = "\t", check.names = FALSE)
  CT <- read.table(file = "Class_taxonomy.txt", sep = "\t", check.names = FALSE)
  OT <- read.table(file = "Order_taxonomy.txt", sep = "\t", check.names = FALSE)
  FT <- read.table(file = "Family_taxonomy.txt", sep = "\t", check.names = FALSE)
  GT <- read.table(file = "Genus_taxonomy.txt", sep = "\t", check.names = FALSE)
  ST <- read.table(file = "Species_taxonomy.txt", sep = "\t", check.names = FALSE)
  
  # Here we are taking the taxonomy tables created above
  # taxa with less than a sum of 0.1 total proportion over all samples are merged
  # into a column of OTHER with the NA column if it is present.
  # NEXT IS TO DYNAMICALLY SET THE 0.1 TO LEAVE 20 TAXA TOTAL WITH OTHER AS REST.
  taxa_list <- list(PT,CT,OT,FT,GT,ST)
  taxa_names <- list("Phylum","Class","Order","Family","Genus","Species")
  for (x in 2:6){
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
    # get top 11 taxa
    keepers <- row.names(as.data.frame(head(sort(colSums(XTtax), decreasing = TRUE), n=11)))
    XT_other <- XTtax[,keepers]
    # sum taxa columns smaller then top 12
    
    if (ncol(XTtax) > 11){
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
    }
    else{XT_other = XTtax}
  
    
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
    colourCount = 12
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    #Sites
    melted$variable = with(melted, reorder(variable,value,mean))
    p2 <- ggplot(data = melted[!is.na(melted$Sites),], ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + facet_grid(~Sites, scales = "free", space='free') + scale_x_discrete(expand = c(0, 0.5)) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) + scale_fill_manual(values = getPalette(colourCount),name = name_label) + ggtitle("Sites")
    filename3 <- paste0(name_label,"Sites_taxonomy.png")
    ggplot2::ggsave(p2, file = filename3, dpi  = 800, width = 8, height = 6, units = "in")
    
    #Times
    melted$variable = with(melted, reorder(variable,value,mean))
    p2 <- ggplot(data = melted[!is.na(melted$Time),], ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + facet_grid(~Time, scales = "free", space='free') + scale_x_discrete(expand = c(0, 0.5)) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) + scale_fill_manual(values = getPalette(colourCount),name = name_label) + ggtitle("Time")
    filename3 <- paste0(name_label,"Times_taxonomy.png")
    ggplot2::ggsave(p2, file = filename3, dpi  = 800, width = 8, height = 6, units = "in")
    
    #GP
    melted$variable = with(melted, reorder(variable,value,mean))
    p2 <- ggplot(data = melted[!is.na(melted$Galleries_Phloem),], ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + facet_grid(~Galleries_Phloem, scales = "free", space='free') + scale_x_discrete(expand = c(0, 0.5)) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) + scale_fill_manual(values = getPalette(colourCount),name = name_label) + ggtitle("Galleries_Phloem")
    filename3 <- paste0(name_label,"Galleries_Phloem_taxonomy.png")
    ggplot2::ggsave(p2, file = filename3, dpi  = 800, width = 8, height = 6, units = "in")
    
    #GA
    melted$variable = with(melted, reorder(variable,value,mean))
    p2 <- ggplot(data = melted[!is.na(melted$Galleries_Adults),], ggplot2::aes(Samples, (value*100), fill = variable)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + facet_grid(~Galleries_Adults, scales = "free", space='free') + scale_x_discrete(expand = c(0, 0.5)) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8)) + scale_fill_manual(values = getPalette(colourCount),name = name_label) + ggtitle("Galleries_Adults")
    filename3 <- paste0(name_label,"Galleries_Adults_taxonomy.png")
    ggplot2::ggsave(p2, file = filename3, dpi  = 800, width = 8, height = 6, units = "in")
    
    
    # heatmap
    #melted$variable = with(melted, reorder(variable,value,mean))
    #p2 <- ggplot(melted, aes(Samples,variable, fill=value)) + geom_tile(color="white", size=0.1) + facet_grid(~Sites, scales = "free_y") + xlab("Day") + ylab(name_label) + theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colours = pal,name = "Proportion") + ggtitle("Sites")
    
  }
}

