#' Diversity Plots
#'
#' This function creates plots that summarize alpha and beta diversity
#' @param brayWmeta table from diversity.R
#' @param newmap table from Create_Tables.R
#' @importFrom magrittr %>%
#'
#' @export

# read in files
Diversity_Plots <- function(brayWmeta,newmap){
  # brayWmeta <- read.table('proportional_diversity_stats.txt')
  # newmap <- read.table("Metadata_common.txt")
  var_explained = (brayWmeta$EV/sum(brayWmeta$EV))*100
  var_explained = format(round(var_explained, 2), nsmall = 2)
  #PCOA <- cmdscale(d2.mes, k = 3, eig = TRUE)
  PCOA$eig <- PCOA$eig
  PCOA$eig <- format(round(PCOA$eig, 2), nsmall = 2)

  Adiv <- function(.data, .column) {
    shannon <- simpson <- chao1 <- NULL
    ShanD <- ggplot2::ggplot(.data, ggplot2::aes(!!dplyr::sym(.column), shannon, colour = !!dplyr::sym(.column))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + ggpubr::stat_compare_means()
    SimD <- ggplot2::ggplot(.data, ggplot2::aes(!!dplyr::sym(.column), simpson, colour = !!dplyr::sym(.column))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Simpson") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + ggpubr::stat_compare_means()
    SimI <- ggplot2::ggplot(.data, ggplot2::aes(!!dplyr::sym(.column), chao1, colour = !!dplyr::sym(.column))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::ylab("Chao1") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + ggpubr::stat_compare_means()
    combinded_plot <- patchwork::wrap_plots(ShanD, SimD, SimI, ncol=1)
    plottitle <- paste0("AlphaDiversity_",.column,".png")
    ggplot2::ggsave(combinded_plot, file=plottitle, dpi=800, height = 12, width = 6, units = "in")
  }
  Adivplots <- names(newmap) %>% purrr::map(~Adiv(.data = brayWmeta, .column = .x))

  Adivnolegend <- function(.data, .column) {
    shannon <- simpson <- chao1 <- NULL
    ShanD <- ggplot2::ggplot(.data, ggplot2::aes(!!dplyr::sym(.column), shannon, colour = !!dplyr::sym(.column))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + ggpubr::stat_compare_means()
    SimD <- ggplot2::ggplot(.data, ggplot2::aes(!!dplyr::sym(.column), simpson, colour = !!dplyr::sym(.column))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Simpson") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + ggpubr::stat_compare_means()
    SimI <- ggplot2::ggplot(.data, ggplot2::aes(!!dplyr::sym(.column), chao1, colour = !!dplyr::sym(.column))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Chao1") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + ggpubr::stat_compare_means()
    combinded_plot <- patchwork::wrap_plots(ShanD, SimD, SimI, ncol=1)
    plottitle <- paste0("AlphaDiversity_",.column,"nolegend.png")
    ggplot2::ggsave(combinded_plot, file=plottitle, dpi=800, height = 12, width = 6, units = "in")
  }
  Adivplots <- names(newmap) %>% purrr::map(~Adivnolegend(.data = brayWmeta, .column = .x))


  # ggplot functions for PCAs
  bdiv <- function(.data, .column) {
    PC1 <- PC2 <- PC3 <- outtab <- NULL
    PC1PC2 <- ggplot2::ggplot(.data, ggplot2::aes(PC1,PC2, colour = !!dplyr::sym(.column))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC2: ",(var_explained[2]), "% variance"))
    PC1PC3 <- ggplot2::ggplot(.data, ggplot2::aes(PC1,PC3, colour = !!dplyr::sym(.column))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC3: ",(var_explained[3]), "% variance"))
    combinded_plot2 <- patchwork::wrap_plots(PC1PC2, PC1PC3, ncol=1)
    plottitle <- paste0("PCA_PC12_PC13_continuous",.column,".png")
    ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
  }
  bdivplots <- names(newmap) %>% purrr::map(~bdiv(.data = brayWmeta, .column = .x))
if("PC1pcoa" %in% colnames(brayWmeta)){
  # ggplot functions for PCoAs
  bdivrare <- function(.data, .column) {
    PC1pcoa <- PC2pcoa <- PC3pcoa <- outtab <- NULL
    PC1PC2 <- ggplot2::ggplot(.data, ggplot2::aes(PC1pcoa,PC2pcoa, colour = !!dplyr::sym(.column))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",PCOA$eig[1], " eigenvalues")) + ggplot2::ylab(paste0("PC2: ",PCOA$eig[2], " eigenvalues"))
    PC1PC3 <- ggplot2::ggplot(.data, ggplot2::aes(PC1pcoa,PC3pcoa, colour = !!dplyr::sym(.column))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")  + ggplot2::xlab(paste0("PC1: ",PCOA$eig[1], " eigenvalues")) + ggplot2::ylab(paste0("PC3: ",PCOA$eig[3], " eigenvalues"))
    combinded_plot2 <- patchwork::wrap_plots(PC1PC2, PC1PC3, ncol=1)
    plottitle <- paste0("Rareified_PCoA_PC12_PC13_continuous",.column,".png")
    ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
  }
  bdivrareplots <- names(newmap) %>% purrr::map(~bdivrare(.data = brayWmeta, .column = .x))
}
  Adiv2 <- function(x) {
    subbed <- brayWmeta[!brayWmeta[,x] == "",]
    brayWmeta <- subbed %>% tidyr::drop_na(x)
    my_comparisons <- unique(brayWmeta[,x])
    out <- tidyr::crossing(my_comparisons,my_comparisons)
    out2 <- subset(out, out[,1] != out[,2])
    makecomp <- out2[as.character(out2$my_comparisons...1) < as.character(out2$my_comparisons...2),]
    comparisons = list(makecomp[1,], makecomp[2,],makecomp[3,],makecomp[4,])
    ShanD <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$shannon, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + ggpubr::stat_compare_means(comparisons = comparisons) + ggpubr::stat_compare_means()
    SimD <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$simpson, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Simpson") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank()) + ggpubr::stat_compare_means(comparisons =comparisons) + ggpubr::stat_compare_means()
    SimI <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$chao1, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::ylab("Chao1") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + ggpubr::stat_compare_means(comparisons = comparisons) + ggpubr::stat_compare_means()
    combinded_plot <- patchwork::wrap_plots(ShanD, SimD, SimI, ncol=1)
    plottitle <- paste0("AlphaDiversity_plots_",x,"_significance.png")
    ggplot2::ggsave(combinded_plot, file=plottitle, dpi=800, height = 12, width = 6, units = "in")
  }

  bdivnolegend <- function(.data, .column) {
    PC1 <- PC2 <- PC3 <- outtab <- NULL
    PC1PC2 <- ggplot2::ggplot(.data, ggplot2::aes(PC1,PC2, colour = !!dplyr::sym(.column))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC2: ",(var_explained[2]), "% variance"))
    PC1PC3 <- ggplot2::ggplot(.data, ggplot2::aes(PC1,PC3, colour = !!dplyr::sym(.column))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC3: ",(var_explained[3]), "% variance"))
    combinded_plot2 <- patchwork::wrap_plots(PC1PC2, PC1PC3, ncol=1)
    plottitle <- paste0("PCA_PC12_PC13_nolegend",.column,".png")
    ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
  }
  bdivplots <- names(newmap) %>% purrr::map(~bdivnolegend(.data = brayWmeta, .column = .x))

# make plots
#lapply(colnames(newmap), Adiv)
#bdiv("Eligibility_group")
#bdiv("V1_ACV2SPIKE_Result")
#lapply(colnames(newmap), bdiv)
}
