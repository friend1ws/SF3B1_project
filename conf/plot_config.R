my_theme <- function() {
  theme_bw(base_family = "Helvetica") %+replace% 
    theme(title = element_text(size = 7),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "grey20", size = 0.5), 
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 7),
          legend.key = element_blank(), 
          legend.key.size = unit(0.25, "cm"),
          legend.margin = margin(0.5, 0.5, 0.5, 0.5),
          legend.text = element_text(size = 6),
          strip.text = element_text(size = 6),
          # strip.background = element_rect(fill = "white", colour = "black", size = 1),
          strip.background = element_blank(), 
          complete = TRUE)
}   