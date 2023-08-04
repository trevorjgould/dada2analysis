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
    ShanD <- ggplot2::ggplot(brayWmeta, ggplot2::aes(x, brayWmeta$shannon, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
    SimD <- ggplot2::ggplot(brayWmeta, ggplot2::aes(x, brayWmeta$simpson, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Simpson") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
    SimI <- ggplot2::ggplot(brayWmeta, ggplot2::aes(x, brayWmeta$chao1, colour = x)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.3)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::ylab("Chao1") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    #combinded_plot <- ShanD / SimD / SimI
    combinded_plot <- gridExtra::grid.arrange(ShanD,SimD,SimI, ncol = 1)
    plottitle <- paste0("AlphaDiversity_plots_",x,".png")
    ggplot2::ggsave(combinded_plot, file=plottitle, dpi=800, height = 12, width = 6, units = "in")
}

# ggplot functions for PCoAs
bdiv <- function(j) {
    PC1PC2 <- ggplot2::ggplot(brayWmeta, ggplot2::aes(brayWmeta$PC1,brayWmeta$PC2, colour = j)) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC2: ",(var_explained[2]), "% variance"))
    PC1PC3 <- ggplot2::ggplot(brayWmeta, ggplot2::aes(brayWmeta$PC1,brayWmeta$PC3, colour = j)) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC3: ",(var_explained[3]), "% variance"))
    #combinded_plot2 <- PC1PC2 / PC1PC3
    combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
    plottitle <- paste0("PCoA_PC12_PC13_continuous",j,".png")
    ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
}

# make plots
lapply(colnames(newmap), Adiv)
lapply(colnames(newmap), bdiv)

##HARDCODED Day Treatment
#ShanD <- ggplot2::ggplot(brayWmeta, ggplot2::aes(Treatment, shannon, colour = as.factor(Day))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.7)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Shannon") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
#SimD <- ggplot2::ggplot(brayWmeta, ggplot2::aes(Treatment, simpson, colour = as.factor(Day))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.7)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA") + ggplot2::ylab("Simpson") + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
#SimI <- ggplot2::ggplot(brayWmeta, ggplot2::aes(Treatment, chao1, colour = as.factor(Day))) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position=ggplot2::position_jitterdodge(),alpha=0.7)+ ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::ylab("Chao1") + ggplot2::xlab("Site") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + scale_color_discrete(name = "Day")
#combinded_plot <- gridExtra::grid.arrange(ShanD,SimD,SimI, ncol = 1)
#plottitle <- paste0("AlphaDiversity_plots_Treatment_Day.png")
#ggplot2::ggsave(combinded_plot, file=plottitle, dpi=800, height = 12, width = 6, units = "in")

# all site combined
PC1PC2 <- ggplot2::ggplot(brayWmeta, ggplot2::aes(PC1,PC2, colour = as.factor(Day), shape = Treatment)) + geom_point(shape=24) + theme_bw() + theme(legend.position = "NA")  + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ylab(paste0("PC2: ",(var_explained[2]), "% variance"))
PC1PC3 <- ggplot2::ggplot(brayWmeta, ggplot2::aes(PC1,PC3, colour = as.factor(Day), shape = Treatment)) + geom_point(shape=24) + theme_bw() + theme(legend.position = "bottom")  + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ylab(paste0("PC3: ",(var_explained[3]), "% variance")) + scale_color_discrete(name = "Day")
#combinded_plot2 <- PC1PC2 / PC1PC3
combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
plottitle <- paste0("PCoA_PC12_PC13_Day_Treatment.png")
ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
# hard coded Site facet species
PC1PC2 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC2, colour = as.factor(brayWmeta$Site))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "NA")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC2: ",(var_explained[2]), "% variance")) + facet_wrap(~brayWmeta$Species)
PC1PC3 <- ggplot2::ggplot(brayWmeta, ggplot2::aes_string(brayWmeta$PC1,brayWmeta$PC3, colour = as.factor(brayWmeta$Site))) + ggplot2::geom_point(size=2) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")  + ggplot2::xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + ggplot2::ylab(paste0("PC3: ",(var_explained[3]), "% variance")) + facet_wrap(~brayWmeta$Species) + labs(color = "Site")
#combinded_plot2 <- PC1PC2 / PC1PC3
combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
plottitle <- paste0("PCoA_PC12_PC13_Site_facet_species.png")
ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 8, width = 6, units = "in")
}
###############
# circle around PCoA PCA points
# ggplot functions for PCoAs
j = "Group"
color4 = c("#BB3717","#54969A","#B47F6D","#0F5372")
color3 = c("#A9502D","#678076","#C7723A")
if(j=="TimePoint") {mycolors = color3}
if(j=="Group") {mycolors = color4}

PC1PC2 <- ggplot(brayWmeta, aes(PC1,PC2, colour = Group, label = ShortID)) + geom_point(size=2) +geom_text_repel() + theme_bw() + theme(legend.position = "NA")  + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + scale_color_manual(values = mycolors) + ylab(paste0("PC2: ",(var_explained[2]), "% variance")) + coord_fixed(ratio = 1)
PC1PC3 <- ggplot(brayWmeta, aes(PC1,PC3, colour = Group, label = ShortID)) + geom_point(size=2) +geom_text_repel() + theme_bw() + theme(legend.position = "bottom")  + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + scale_color_manual(values = mycolors) + ylab(paste0("PC3: ",(var_explained[3]), "% variance")) + coord_fixed(ratio = 1)
#combinded_plot2 <- PC1PC2 / PC1PC3
combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
plottitle <- paste0("PCoA_PC12_PC13_continuous",j,"facet_timepoint.png")
#plottitle <- paste0("PCoA_PC12_PC13_continuous",j,".png")
ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 6, width = 12, units = "in")

# cage and Group labelled
color4 = c("#BB3717","#54969A","#B47F6D","#0F5372")
PC1PC2 <- ggplot(brayWmeta, aes(PC1,PC2, color = as.factor(Cage), label = ShortID)) + geom_point(size=2, aes(shape = Group)) +geom_text_repel() + theme_bw() + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + scale_color_manual(values = color4,name = "Cage") + ylab(paste0("PC2: ",(var_explained[2]), "% variance")) + coord_fixed(xlim = c(-15,15),ylim = c(-15,15)) + stat_ellipse(aes(group=as.factor(Cage), type = "norm"))
PC1PC3 <- ggplot(brayWmeta, aes(PC1,PC3, color = as.factor(Cage), label = ShortID)) + geom_point(size=2, aes(shape = Group)) +geom_text_repel() + theme_bw() + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + scale_color_manual(values = color4,name = "Cage") + ylab(paste0("PC3: ",(var_explained[3]), "% variance")) + coord_fixed(xlim = c(-15,15),ylim = c(-15,15)) + stat_ellipse(aes(group=as.factor(Cage), type = "norm"))
#combinded_plot2 <- PC1PC2 / PC1PC3
combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
plottitle <- paste0("PCA_PC12_PC13_continuous_Group_Cage_facet_timepoint_labelled.png")
ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 12, width = 6, units = "in")   

# cage and Group NOT labelled
color4 = c("#BB3717","#54969A","#B47F6D","#0F5372")
PC1PC2 <- ggplot(brayWmeta, aes(PC1,PC2, color = as.factor(Cage))) + geom_point(size=2, aes(shape = Group)) + theme_bw() + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + scale_color_manual(values = color4,name = "Cage") + ylab(paste0("PC2: ",(var_explained[2]), "% variance")) + coord_fixed(xlim = c(-15,15),ylim = c(-15,15)) + stat_ellipse(aes(group=as.factor(Cage), type = "norm"))
PC1PC3 <- ggplot(brayWmeta, aes(PC1,PC3, color = as.factor(Cage))) + geom_point(size=2, aes(shape = Group)) + theme_bw() + xlab(paste0("PC1: ",(var_explained[1]), "% variance")) + scale_color_manual(values = color4,name = "Cage") + ylab(paste0("PC3: ",(var_explained[3]), "% variance")) + coord_fixed(xlim = c(-15,15),ylim = c(-15,15)) + stat_ellipse(aes(group=as.factor(Cage), type = "norm"))
#combinded_plot2 <- PC1PC2 / PC1PC3
combinded_plot2 <- gridExtra::grid.arrange(PC1PC2, PC1PC3, ncol = 1)
plottitle <- paste0("PCA_PC12_PC13_continuous_Group_Cage_facet_timepoint_NOT_labelled.png")
ggplot2::ggsave(combinded_plot2, file=plottitle, dpi=800, height = 12, width = 6, units = "in") 


