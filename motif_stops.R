library(ggplot2)
library(ggarrange)

stop_count_plot <- function(file_path, main=NULL) {
  library(ggplot2)
  file <- read.csv(file_path, sep=",", head=T)
  simulants <- file[file$id != 'real',]
  real <- file[file$id == 'real',]
  
  plot <- ggplot(data=simulants, aes(simulants$stop_count)) +
    geom_histogram(breaks=seq(min(file$stop_count), max(file$stop_count), by = 1), col="black", fill="RoyalBlue", alpha = 1) +
    labs(x="Stop codon count", y="Frequency") +
    geom_vline(xintercept=file$stop_count[file$id == "real"], lty=2, col="red") + 
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
  
  greater_than <- (nrow(simulants[simulants$stop_count <= real$stop_count,]) + 1) / (nrow(simulants) + 1)
  print(greater_than)
  
  return(plot)
}

get_z <- function(data){
  real = data[data$simulation == "real",]
  sims = data[data$simulation != "real",]
  z <- (real$count - mean(sims$count)) / sd(sims$count)
  return(z)
}

p_from_z <- function(z) {
  return(2*pnorm(-abs(z)))
}


# motif stop counts
eses <- stop_count_plot("results/motif_stops_simulation/CaceresHurstESEs_INT3/stop_counts_10000.txt", "INT3 ESEs")
ises <- stop_count_plot("results/motif_stops_simulation/ises_wang_2012/stop_counts_10000.txt", "Wang 2012 ISEs")
plot <- ggarrange(eses, ises, ncol=2, labels=c("A", "B"))
ggsave('results/graphs/motif_stop_counts.pdf', width=10, height=5)
ggsave('results/graphs/motif_stop_counts.eps', width=10, height=5)
