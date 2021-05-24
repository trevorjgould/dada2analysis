#' Diversity Plots
#'
#' This function creates plots that summarize alpha and beta diversity
#' @param brayWmeta table from diversity.R
#' @param newmap table from Create_Tables.R
#'
#' @export
#' @importFrom gridExtra grid.arrange

# read in files
Diversity_Plots <- function(brayWmeta,newmap){
# brayWmeta <- read.table('proportional_diversity_stats.txt')
# newmap <- read.table("Metadata_common.txt")
var_explained = (brayWmeta$EV/sum(brayWmeta$EV))*100
var_explained = format(round(var_explained, 2), nsmall = 2)

Adiv <- function(x) {
    ShanD <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$shannon, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
    SimD <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$simpson, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Simpson") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
    SimI <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(x, brayWmeta$invsimpson, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::ylab("Inverse_Simpson") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    #combinded_plot <- ShanD / SimD / SimI
    combinded_plot <- gridExtra::grid.arrange(ShanD,SimD,SimI, ncol = 1)
    plottitle <- paste0("AlphaDiversity_plots_",x,".png")
    ggplot2::ggsave(combinded_plot, file=plottitle, dpi=800, height = 12, width = 6, units = "in")
}

# ggplot functions for PCoAs
bdiv <- function(j) {
    PC1PC2 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC2, colour = j)) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC2: ",(var_explained[2]), "% variance"))
    PC1PC3 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC3, colour = j)) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC3: ",(var_explained[3]), "% variance"))
    #combinded_plot2 <- PC1PC2 / PC1PC3
    combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
    plottitle <- paste0("PCoA_PC12_PC13_continuous",j,".png")
    ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
}

# make plots
lapply(colnames(newmap), Adiv)
lapply(colnames(newmap), bdiv)

#HARDCODED SITE PLOT
ShanD <- ggplot2::ggplot(brayWmeta, ggplot2::aes(as.factor(Site), shannon, colour = as.factor(Site))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
SimD <- ggplot2::ggplot(brayWmeta, ggplot2::aes(as.factor(Site), simpson, colour = as.factor(Site))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
SimI <- ggplot2::ggplot(brayWmeta, ggplot2::aes(as.factor(Site), invsimpson, colour = as.factor(Site))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Inverse_Simpson") + ggplot2::xlab("Site") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
#combinded_plot <- ShanD / SimD / SimI
combinded_plot <- gridExtra::grid.arrange(ShanD,SimD,SimI, ncol = 1)
plottitle <- paste0("AlphaDiversity_plots_Site.png")
ggplot2::ggsave(combinded_plot, file=plottitle, dpi=800, height = 12, width = 6, units = "in")
# all site combined
PC1PC2 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC2, colour = as.factor(brayWmeta$Site))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC2: ",(var_explained[2]), "% variance"))
PC1PC3 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC3, colour = as.factor(brayWmeta$Site))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC3: ",(var_explained[3]), "% variance"))
#combinded_plot2 <- PC1PC2 / PC1PC3
combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
plottitle <- paste0("PCoA_PC12_PC13_Site.png")
ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
# hard coded Site facet species
PC1PC2 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC2, colour = as.factor(brayWmeta$Site))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC2: ",(var_explained[2]), "% variance")) + facet_wrap(~brayWmeta$Species)
PC1PC3 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC3, colour = as.factor(brayWmeta$Site))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC3: ",(var_explained[3]), "% variance")) + facet_wrap(~brayWmeta$Species) + labs(color = "Site")
#combinded_plot2 <- PC1PC2 / PC1PC3
combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
plottitle <- paste0("PCoA_PC12_PC13_Site_facet_species.png")
ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
