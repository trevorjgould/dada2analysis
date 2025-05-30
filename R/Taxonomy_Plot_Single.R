#' Create_Taxonomy_Plot
#'
#' This function creates plots that summarize prevalent taxa
#'
#' table of metadata.
#' @param metadata table
#' @param taxtable table
#' @param t integer
#'
#' @export

# here we are looking to create a function that runs a single level table instead of reading in each table that was created

Create_Taxonomy_Plot <- function(taxtable,t,metadata){
	value <- variable <- Samples <- meta <- taxonomy <- NULL
	# get taxa output name
    name_label <- deparse1(substitute(taxtable))
    name_label <- gsub(".txt","",name_label)
    name_label <- gsub(".*\\$", "", name_label)
    # get taxa table
    xt=as.data.frame(taxtable)
    xt=xt[(rowSums(xt)>0),]

    # fix rownames . to -
    # row.names(metadata) <- gsub("-", "\\.", row.names(metadata))

    #find common names
    common <- intersect(rownames(metadata),rownames(xt))
    # get just the overlapping samples
    metadata <- metadata[common,, drop = FALSE]
    xt <- xt[common,, drop = FALSE]


    XTtax = sweep(xt, 1, rowSums(xt),'/')
    # KEEP TOP TAXA AND SORT REST INTO OTHER CATEGORY AND PLOT
    if (ncol(XTtax) <= t){
    XT_other = XTtax
    }
 
    # if number of columns is = or less than summary asked for (t) keep table as is. 
    if (ncol(XTtax) <= t){
 		 XT_other = XTtax
 	}
    # get taxa that are present > t
    if (ncol(XTtax) > t){
      topt <- rownames(as.data.frame(utils::head(sort(colSums(XTtax), decreasing = TRUE), n=t)))
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

    both <- cbind(XT_other,metadata)
    both$Samples <- row.names(metadata)
    melted <- reshape2::melt(both, id.vars = c(colnames(metadata),"Samples"))
    melted$variable <- gsub("[a-z]__", "", melted$variable)
    melted$variable <- gsub("\\.", " ", melted$variable)

    melted$variable = with(melted, reorder(variable,value,mean))
    melted$SampleID <- rownames(melted)

    alltaxa <- unique(melted$variable)
    alltaxa2 <- alltaxa[alltaxa != "other"]
    melted$taxonomy <- factor(melted$variable, levels = c("other",paste0(alltaxa2)))
    # filenames
    filename <- paste0(name_label,"_taxonomy_top",t,"_Type_boxplot.png")
    filename2 <- paste0(name_label,"_taxonomy_top_",t,".txt")
    titlename <- paste0("Top_",t,"_",name_label)
    # colors
    c25 <- c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
    AsteroidCityall = c("#0A9F9D", "#CEB175", "#E54E21", "#6C8645", "#C18748","#C52E19", "#AC9765", "#54D8B1", "#b67c3b", "#175149", "#AF4E24", "#FBA72A", "#D3D4D8", "#CB7A5C", "#5785C1")
	GrandBudapest1 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4","#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
	#getPalette = RColorBrewer::colorRampPalette(brewer.pal(t, "Set1"))

	p <- ggplot2::ggplot(data = melted[!is.na(melted$SampleID),], ggplot2::aes(Samples, (value*100), fill = taxonomy)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + ggplot2::scale_x_discrete(expand = c(0, 0.5)) + ggplot2::xlab("Samples") + theme_Publication() + ggplot2::scale_fill_manual(values = c25) + ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2,title=name_label))
    ggplot2::ggsave(p, file = filename, dpi  = 800, width = 14, height = 8, units = "in")
    utils::write.table(both, file = filename2, sep = "\t", quote = FALSE)
    return(p)
}

Create_Taxonomy_Plot_Facet <- function(taxtable,t,metadata){
    value <- variable <- Samples <- meta <- taxonomy <- NULL
    # get taxa output name
    facets <- colnames(metadata)
    name_label <- deparse1(substitute(taxtable))
    name_label <- gsub(".txt","",name_label)
    name_label <- gsub(".*\\$", "", name_label)
    # get taxa table
    xt=as.data.frame(taxtable)
    xt=xt[(rowSums(xt)>0),]
    
    # fix rownames . to -
    # row.names(metadata) <- gsub("-", "\\.", row.names(metadata))
    
    #find common names
    common <- intersect(rownames(metadata),rownames(xt))
    # get just the overlapping samples
    metadata <- metadata[common,, drop = FALSE]
    xt <- xt[common,, drop = FALSE]
    
    
    XTtax = sweep(xt, 1, rowSums(xt),'/')
    # KEEP TOP TAXA AND SORT REST INTO OTHER CATEGORY AND PLOT
    if (ncol(XTtax) <= t){
    XT_other = XTtax
    }
 
    # get taxa that are present > t
    if (ncol(XTtax) > t){
        topt <- rownames(as.data.frame(utils::head(sort(colSums(XTtax), decreasing = TRUE), n=t)))
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
    
    both <- cbind(XT_other,metadata)
    both$Samples <- row.names(metadata)
    melted <- reshape2::melt(both, id.vars = c(colnames(metadata),"Samples"))
    melted$variable <- gsub("[a-z]__", "", melted$variable)
    melted$variable <- gsub("\\.", " ", melted$variable)
    
    melted$variable = with(melted, reorder(variable,value,mean))
    melted$SampleID <- rownames(melted)
    
    alltaxa <- unique(melted$variable)
    alltaxa2 <- alltaxa[alltaxa != "other"]
    melted$taxonomy <- factor(melted$variable, levels = c("other",paste0(alltaxa2)))
    # filenames

    # colors
    c25 <- c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
    AsteroidCityall = c("#0A9F9D", "#CEB175", "#E54E21", "#6C8645", "#C18748","#C52E19", "#AC9765", "#54D8B1", "#b67c3b", "#175149", "#AF4E24", "#FBA72A", "#D3D4D8", "#CB7A5C", "#5785C1")
    GrandBudapest1 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4","#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
    #getPalette = RColorBrewer::colorRampPalette(brewer.pal(t, "Set1"))
facetPlot <- function(i){
		filename <- paste0(name_label,"_taxonomy_top",t,"_in_",i,"_Type_boxplot.png")
    	filename2 <- paste0(name_label,"_taxonomy_top_",t,"_in_",i,".txt")
   		titlename <- paste0("Top_",t,"_in_",i,"_",name_label)
        p <- ggplot2::ggplot(data = melted[!is.na(melted$SampleID),], ggplot2::aes(Samples, (value*100), fill = taxonomy)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + ggplot2::scale_x_discrete(expand = c(0, 0.5)) + ggplot2::xlab("Samples") + ggplot2::facet_wrap(~get(i), scales = "free") + theme_Publication() + ggplot2::scale_fill_manual(values = c25) + ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2,title=name_label))
        ggplot2::ggsave(p, file = filename, dpi  = 800, width = 14, height = 8, units = "in")
        utils::write.table(both, file = filename2, sep = "\t", quote = FALSE)
        return(p)
    }
    
    plotlist <- lapply(facets,facetPlot)
}

Create_Taxonomy_Plot_Facet_single <- function(taxtable,t,metadata,var1){
    value <- variable <- Samples <- meta <- taxonomy <- NULL
    # get taxa output name
    facets <- var1
    name_label <- deparse1(substitute(taxtable))
    name_label <- gsub(".txt","",name_label)
    name_label <- gsub(".*\\$", "", name_label)
    # get taxa table
    xt=as.data.frame(taxtable)
    xt=xt[(rowSums(xt)>0),]
    
    # fix rownames . to -
    # row.names(metadata) <- gsub("-", "\\.", row.names(metadata))
    
    #find common names
    common <- intersect(rownames(metadata),rownames(xt))
    # get just the overlapping samples
    metadata <- metadata[common,, drop = FALSE]
    xt <- xt[common,, drop = FALSE]
    
    
    XTtax = sweep(xt, 1, rowSums(xt),'/')
    # KEEP TOP TAXA AND SORT REST INTO OTHER CATEGORY AND PLOT
    if (ncol(XTtax) <= t){
    XT_other = XTtax
    }
 
    # get taxa that are present > t
    if (ncol(XTtax) > t){
        topt <- rownames(as.data.frame(utils::head(sort(colSums(XTtax), decreasing = TRUE), n=t)))
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
    
    both <- cbind(XT_other,metadata)
    both$Samples <- row.names(metadata)
    melted <- reshape2::melt(both, id.vars = c(colnames(metadata),"Samples"))
    melted$variable <- gsub("[a-z]__", "", melted$variable)
    melted$variable <- gsub("\\.", " ", melted$variable)
    
    melted$variable = with(melted, reorder(variable,value,mean))
    melted$SampleID <- rownames(melted)
    
    alltaxa <- unique(melted$variable)
    alltaxa2 <- alltaxa[alltaxa != "other"]
    melted$taxonomy <- factor(melted$variable, levels = c("other",paste0(alltaxa2)))
    # filenames

    # colors
    c25 <- c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
    AsteroidCityall = c("#0A9F9D", "#CEB175", "#E54E21", "#6C8645", "#C18748","#C52E19", "#AC9765", "#54D8B1", "#b67c3b", "#175149", "#AF4E24", "#FBA72A", "#D3D4D8", "#CB7A5C", "#5785C1")
    GrandBudapest1 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4","#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
    #getPalette = RColorBrewer::colorRampPalette(brewer.pal(t, "Set1"))
facetPlot <- function(i){
		filename <- paste0(name_label,"_taxonomy_top",t,"_in_",i,"_Type_boxplot.png")
    	filename2 <- paste0(name_label,"_taxonomy_top_",t,"_in_",i,".txt")
   		titlename <- paste0("Top_",t,"_in_",i,"_",name_label)
        p <- ggplot2::ggplot(data = melted[!is.na(melted$SampleID),], ggplot2::aes(Samples, (value*100), fill = taxonomy)) + ggplot2::geom_bar(stat='identity')+ ggplot2::ylab("Percent") + ggplot2::scale_x_discrete(expand = c(0, 0.5)) + ggplot2::xlab("Samples") + ggplot2::facet_wrap(~get(i), scales = "free") + theme_Publication() + ggplot2::scale_fill_manual(values = c25) + ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2,title=name_label))
        ggplot2::ggsave(p, file = filename, dpi  = 800, width = 14, height = 8, units = "in")
        utils::write.table(both, file = filename2, sep = "\t", quote = FALSE)
        return(p)
    }
    
    facetPlot(facets)
}
