library("ggplot2")

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0) {
  stop("At least one argument must be supplied -number of simulations", call.=FALSE)
} else {
   simulations <- args[1]
}

# set up blank list

if (args[2] == "all_sets") {
  sets <- c(
    # "RBP_motifs/RBP_motifs_stop_counts_pos_nd",
    # "RBP_motifs/RBP_motifs_stop_counts_neg_nd",
    "ESR/ESR_stop_counts",
    "INT3/INT3_stop_counts",
    "Ke400_ESEs/Ke400_ESEs_stop_counts",
    "PESE/PESE_stop_counts",
    "RBP_motifs/RBP_motifs_stop_counts",
    "RESCUE/RESCUE_stop_counts"
  )
} else {
  sets <- c()
  for (arg in 2:length(args)){
    sets <- c(sets, args[arg])
  }
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


get_plot <- function(file_path, out_path, main=NULL, title_left=NULL) {
  file <- read.csv(file_path, sep=",", head=T)
  simulants <- file[file$id != 'real',]
  real <- file[file$id == 'real',]
  
  plot <- ggplot(data=simulants, aes(simulants$stop_count)) +
    geom_histogram(breaks=seq(min(file$stop_count), max(file$stop_count), by = 1), col="black", fill="RoyalBlue", alpha = .8) +
    # labs(title=paste(sprintf(main), " (", nrow(simulants), " simulations)", sep="")) +
    labs(x="Stop codon count in all reading frames", y="Frequency") +
    geom_vline(xintercept=file$stop_count[file$id == "real"], lty=2, col="red") + 
    plot_theme()
  
  plot <- title_styling(plot, main, title_left)
  # out_path = paste("ese_stops_simulations_output/", set, "/", set, "_", simulations, "_hist.head(id", sep="")
  # ggsave(out_path, plot=plot)
  print(plot)
  # return(plot)
}

get_plot("results/motif_stops_simulation/CaceresHurstESEs_INT3/stop_counts_10000.txt", "results/motif_stops_simulation/CaceresHurstESEs_INT3/plot.pdf", "ESEs INT3", "A")


run_plots <- function(simulations) {
  for (set in sets) {
    print(set)
    file_path <- paste("ese_stops_simulations_output/", set, "/", set, "_stop_counts_" , simulations, ".csv", sep="")
    out_path <- paste("ese_stops_simulations_output/", set, "/", set, "_hist.jpg", sep="")
    get_plot(file_path, out_path)
  }
}

run_plots(simulations)

prepare_data <- function(set) {
  if (set == "grouped") {
    sample_file <- read.table('results/ese_stops_simulations_output/RBP_motifs_grouped/motif_samples.txt', head=T)
  }
}

# sample_file <- read.table('results/ese_stops_simulations_output/RBP_motifs_grouped/motif_samples.txt')
# colnames(sample_file) <- c("sample")
# 
# for (set in sample_file$sample) {
#   file_path = paste("results/ese_stops_simulations_output/RBP_motifs_grouped/", set, "/", set, "_stop_counts_1000.csv", sep="")
#   out_path = paste("results/ese_stops_simulations_output/RBP_motifs_grouped/plots/", set, "_hist.jpg", sep="")
#   get_plot(file_path, out_path)
# }


get_p <- function(file_path, direction="greater") {
  
  file <- read.table(file_path, head=T, sep=",")
  real_count = file[file$id == "real",]
  simulants = file[file$id != "real",]
  if (direction == "less") {
    p <- (nrow(simulants[simulants$stop_count < real_count$stop_count,] + 1)) / nrow(simulants) 
  } else {
    p <- (nrow(simulants[simulants$stop_count > real_count$stop_count,] + 1)) / nrow(simulants)
  }
  return(p)
}

get_plot("results/motif_stops_simulation/CaceresHurstESEs_INT3/stop_counts_10000.txt", "results/motif_stops_simulation/CaceresHurstESEs_INT3/plot.pdf", "A", title_left = T)
get_plot("results/motif_stops_simulation/ess_fas-hex3_wang_2004/stop_counts_10000.txt", "results/motif_stops_simulation/ess_fas-hex3_wang_2004/plot.pdf", "B", title_left=T)
get_plot("results/motif_stops_simulation/ises_wang_2012/stop_counts_10000.txt", "results/motif_stops_simulation/ises_wang_2012/plot.pdf", "B", title_left=T)

get_p("results/motif_stops_simulation/CaceresHurstESEs_INT3/stop_counts_10000.txt", direction="less")
get_p("results/motif_stops_simulation/ises_wang_2012/stop_counts_10000.txt", direction="less")
get_p("results/motif_stops_simulation/ess_fas-hex3_wang_2004/stop_counts_10000.txt", direction="less")

# get_p("results/motif_stops_simulation/ises_wang_2012/stop_counts_10000.txt", direction="less")
# get_plot("results/motif_stops_simulation/ises_wang_2012/stop_counts_10000.txt", "results/motif_stops_simulation/ises_wang_2012/plot.pdf", "ISEs Wang 2012")
# 
# get_p("results/motif_stops_simulation/ess_fas-hex3_wang_2004/stop_counts_10000.txt")
# get_plot("results/motif_stops_simulation/ess_fas-hex3_wang_2004/stop_counts_10000.txt", "results/motif_stops_simulation/ess_fas-hex3_wang_2004/plot.pdf", "ESSs Wang 2004")
# 
# get_p("results/motif_stops_simulation/CaceresHurstESEs_INT3/stop_counts_10000.txt", direction="less")
# get_plot("results/motif_stops_simulation/CaceresHurstESEs_INT3/stop_counts_10000.txt", "results/motif_stops_simulation/CaceresHurstESEs_INT3/plot.pdf", "ESEs INT3")
# 
# get_p("results/motif_stops_simulation/filtered_RBP_motifs_nonCDS/stop_counts_10000.txt", direction="less")
# get_plot("results/motif_stops_simulation/filtered_RBP_motifs_nonCDS/stop_counts_10000.txt", "results/motif_stops_simulation/filtered_RBP_motifs_nonCDS/plot.pdf", "RBP motifs non-CDS")
# 
# get_p("results/motif_stops_simulation/ess_wang_2007/stop_counts_10000.txt", direction="less")
# get_plot("results/motif_stops_simulation/ess_wang_2007/stop_counts_10000.txt", "results/motif_stops_simulation/ess_wang_2007/plot.pdf", "ESS Wang 2007")
# 
# get_plot("results/motif_stops_simulation/CaceresHurstESEs_INT3/stop_counts_1000.txt", "results/motif_stops_simulation/ess_wang_2007/plot.pdf", "ESS Wang 2007")
# get_plot("results/motif_stops_simulation/ess_fas-hex3_wang_2004/stop_counts_1000.txt", "results/motif_stops_simulation/ess_wang_2007/plot.pdf", "ESS Wang 2007")
# get_plot("results/motif_stops_simulation/ises_wang_2012/stop_counts_1000.txt", "results/motif_stops_simulation/ess_wang_2007/plot.pdf", "ESS Wang 2007")
# 

# P value for INT3
file_path <- paste("ese_stops_simulations_output/", "INT3", "/", "INT3", "_stop_counts_" , 10000, ".csv", sep="")
file <- read.csv(file_path, head=T)

simulants = file[file$id != "real",]
real = file[file$id == "real",]


(nrow(simulants[simulants$stop_count < real$stop_count,]) + 1) / (nrow(simulants) + 1)
