theme_plt <- function(base_size = 10, base_family = "sans") {
  theme_minimal() %+replace%
    theme(axis.ticks = element_line(colour = "black", size = 0.2),
          panel.border = element_rect(fill = NA),
          panel.grid.major = element_line(colour = "gray80", size = 0.25, linetype = "dashed"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 8,
                                   family = "sans",
                                   color = "gray15"))
}


theme_plt2 <- function(base_size = 9, base_family = "sans") {
  theme_bw(
    base_family = base_family,
    base_size = base_size
  ) +
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(size = 0.4),
      #axis.ticks = element_line(size = 0.3),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = rel(0.9)),
      strip.placement = "outside",
      # strip.background = element_rect(fill = "gray95", color = NA),
      panel.spacing = unit(1.5, "lines"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.text = element_text(size = 9),
      legend.text.align = 0,
      legend.key = element_blank()
    )
}