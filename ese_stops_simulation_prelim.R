library("ggplot2")
setwd(getwd())

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

run_plots <- function(simulations) {
  for (set in sets) {
    print(set)
    file_path <- paste("ese_stops_simulations_output/", set, "/", set, "_stop_counts_" , simulations, ".csv", sep="")
    file <- read.csv(file_path, sep=",", head=T)
    simulants <- file[file$id != 'real',]
    real <- file[file$id == 'real',]

    plot <- ggplot(data=simulants, aes(simulants$stop_count)) +
      geom_histogram(breaks=seq(min(file$stop_count), max(file$stop_count), by = 2), col="black", fill="white", alpha = .8) +
      labs(title=paste(sprintf(set), " (", nrow(simulants), " simulations)", sep="")) +
      labs(x="Stop codon count in all frames", y="Count") +
      geom_vline(xintercept=file$stop_count[file$id == "real"], lty=2, col="red")

    out_path = paste("ese_stops_simulations_output/", set, "/", set, "_", simulations, "_hist.head(id", sep="")
    ggsave(out_path, plot=plot)
  }
}
run_plots(simulations)
