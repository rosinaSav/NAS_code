library(ggplot2)
library(dplyr)

title_styling <- function(plot, title, title_left) {
  if(!is.null(title)) {
    plot <- plot +
      labs(title=title) +
      theme(
        plot.margin=unit(c(10,5,5,5),"mm")
      )
  }
  if(!is.null(title_left)) {
    plot <- plot +
      theme(
        plot.title = element_text(hjust = -0.075)
      )
  }
  return(plot)
}



plot_theme <- function(font_size=12) {
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    text = element_text(),
    panel.background = element_blank(),
    plot.background = element_rect(colour = NA),
    # panel.border = element_rect(colour = NA),
    axis.title = element_text(face = "bold",size = rel(1.2)),
    axis.title.y = element_text(angle=90,vjust =2),
    axis.title.x = element_text(vjust = -0.2),
    axis.text = element_text(size=font_size),
    # axis.text.x = element_text(size=font_size),
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    panel.grid.major = element_line(colour="#f0f0f0"),
    panel.grid.minor = element_blank(),
    # legend.key = element_rect(colour = NA),
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    # legend.key.size= unit(0.2, "cm"),
    # legend.title = element_text(face="italic"),
    plot.margin=unit(c(10,5,5,5),"mm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
  )
}


get_n_t_value <- function(neg_control, feature1, feature2, threshold) {
  n_t = rep(NA, length(neg_control))
  diffs <- data.frame(true = double(), min = double(), max=double())
  for (index in 1:length(neg_control)) {
    # get the true difference
    true_diff = NAS_data[index, feature1] - NAS_data[index, feature2]
    sim_diffs = neg_control[[index]][feature1] - neg_control[[index]][feature2]
    sim_diffs = sim_diffs[!is.na(sim_diffs)]
    diffs <- rbind(diffs, data.frame(true_diff, min(sim_diffs), max(sim_diffs)))
  }
  return(diffs)
}

get_n_t_plot_line =  function(feature1, feature2, neg_control, NAS_data, xlab, ylab, swap = FALSE, big_only = FALSE, return_p = TRUE, return_plot = FALSE, title=NULL, title_left=NULL) {
  threshold = 0
  y_pos = 60
  if (big_only != FALSE) {
    threshold = big_only
    y_pos = 10
  }
  if (swap == TRUE) {
    temp = feature1
    feature1 = feature2
    feature2 = temp
  }
  
  n_t <- get_n_t_value(neg_control, feature1, feature2, threshold)
  colnames(n_t) <- c("true_diff", "min", "max")
  n_t$ID <- seq.int(nrow(n_t))
  
  head(n_t)
  
  plot <- ggplot() +
    geom_line(aes(x = n_t$ID, y = n_t$true_diff), col="RoyalBlue") +
    geom_ribbon(aes(x=n_t$ID, ymin = n_t$min, ymax = n_t$max)) +
    geom_hline(yintercept = 0.025, lty=2, col="#aaaaaa") +
    geom_hline(yintercept = -0.025, lty=2, col="#aaaaaa")
    # coord_cartesian(ylim=c(-0.5,0.5))
    
    # scale_y_continuous(limits = c(-20, 100), labels = seq(-20, 100, 10), breaks = seq(-20, 100, 10)) +
    # labs(x = xlab, y=ylab) +
    # theme(
    #   axis.text.x = element_blank(),
    #   axis.ticks.x=element_blank(),
    #   panel.background = element_rect( fill = "#f5f5f5"),
    #   panel.grid.major = element_line(colour="#ffffff"),
    #   panel.grid.minor = element_line(colour="#ffffff"),
    # )

  plot  
}
get_n_t_plot_line("norm_count_no_PTC_incl", "norm_count_het_PTC_incl", neg_control, NAS_data, xlab="", ylab="PSIdiff", title="B", title_left =T)

head(NAS_data)



