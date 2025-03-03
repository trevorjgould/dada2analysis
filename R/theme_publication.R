# Source: https://github.com/koundy/ggplot_theme_Publication/blob/master/ggplot_theme_Publication-2.R

theme_Publication <- function(base_size=16, base_family="Palatino Linotype") {
    
    (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
        + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                size = ggplot2::rel(1.2), hjust = 0.5, margin = ggplot2::margin(0,0,20,0)),
                text = ggplot2::element_text(color="black"),
                panel.background = ggplot2::element_rect(colour = NA),
                plot.background = ggplot2::element_rect(colour = NA),
                panel.border = ggplot2::element_rect(colour = NA),
                axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                axis.title.x = ggplot2::element_text(vjust = -0.2),
                axis.text = ggplot2::element_text(), 
                axis.text.x=ggplot2::element_blank(),
                axis.ticks.x=ggplot2::element_blank(),
                axis.line.x = ggplot2::element_line(colour="black"),
                axis.line.y = ggplot2::element_line(colour="black"),
                axis.ticks = ggplot2::element_line(),
                panel.grid.major = ggplot2::element_line(colour="#f0f0f0"),
                panel.grid.minor = ggplot2::element_blank(),
                legend.key = ggplot2::element_rect(colour = NA),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.box = "vetical",
                legend.key.size= ggplot2::unit(0.5, "cm"),
                #legend.margin = ggplot2::unit(0, "cm"),
                legend.title = ggplot2::element_blank(),
                legend.text=ggplot2::element_text(color="black"),
                plot.margin=ggplot2::unit(c(2,2,2,2),"mm"),
                strip.background=ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = ggplot2::element_text(face="bold")
        ))
    
}

#scale_fill_Publication <- function(...){
#    scales::discrete_scale("fill","Publication",scales::manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
#    
#}
#
#scale_colour_Publication <- function(...){
#    scales::discrete_scale("colour","Publication",scales::manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
#    
#}

theme_Publication_title <- function(base_size=16, base_family="Palatino Linotype") {
    
    (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
        + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                size = ggplot2::rel(1.2), hjust = 0.5, margin = ggplot2::margin(0,0,20,0)),
                text = ggplot2::element_text(color="black"),
                panel.background = ggplot2::element_rect(colour = NA),
                plot.background = ggplot2::element_rect(colour = NA),
                panel.border = ggplot2::element_rect(colour = NA),
                axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                axis.title.x = ggplot2::element_text(vjust = -0.2),
                axis.text = ggplot2::element_text(), 
                axis.text.x=ggplot2::element_blank(),
                axis.ticks.x=ggplot2::element_blank(),
                axis.line.x = ggplot2::element_line(colour="black"),
                axis.line.y = ggplot2::element_line(colour="black"),
                axis.ticks = ggplot2::element_line(),
                panel.grid.major = ggplot2::element_line(colour="#f0f0f0"),
                panel.grid.minor = ggplot2::element_blank(),
                legend.key = ggplot2::element_rect(colour = NA),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.box = "vetical",
                legend.key.size= ggplot2::unit(0.5, "cm"),
                #legend.margin = ggplot2::unit(0, "cm"),
                legend.text=ggplot2::element_text(color="black"),
                plot.margin=ggplot2::unit(c(2,2,2,2),"mm"),
                strip.background=ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = ggplot2::element_text(face="bold")
        ))
    
}
