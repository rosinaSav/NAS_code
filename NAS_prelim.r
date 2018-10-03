library(grid)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

compare_the_two_measures = function(NAS_data, title) {
  cor_test = cor.test(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, NAS_data$norm_count_no_PTC - NAS_data$norm_count_het_PTC, method = "spearman")
  plot(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, NAS_data$norm_count_no_PTC - NAS_data$norm_count_het_PTC, pch = 19, col = "RoyalBlue", xlab = "PSI(-/-) - PSI(-/+)", ylab = "RPMskip(-/-) - RPMskip(-/+)", main = title)
  abline(h = 0)
  abline(v = 0)
  text(x = 10, y = 1, labels = paste("rho ~ ", round(cor_test$estimate, 3), "\np ~ ", round(cor_test$p.value, 3), sep = ""), adj = 0, cex = 1.3)
}

expression_analysis = function(overlap, all_exons, expression) {
  overlap = unlist(lapply(strsplit(overlap, ".", fixed = TRUE), `[[`, 1))
  all_exons = unlist(lapply(strsplit(all_exons, ".", fixed = TRUE), `[[`, 1))
  all_exons = setdiff(all_exons, overlap)
  overlap = unique(overlap)
  for (column in colnames(expression[2:length(colnames(expression))])) {
    expression_analysis_core(overlap, all_exons, column, expression)
  }
}

expression_analysis_core = function(overlap, all_exons, column, expression) {
  print("***")
  print(column)
  overlap_data = expression[expression$gene %in% overlap, column]
  all_data = expression[expression$gene %in% all_exons, column]
  print("Big effects median:")
  print(median(overlap_data, na.rm = TRUE))
  print("Median for everything else:")
  print(median(all_data, na.rm = TRUE))
  print("p for difference:")
  print(wilcox.test(overlap_data, all_data, alt = "t")$p.value)
}


# get_n_t_value(neg_control, "PSI_no_PTC", "PSI_het_PTC", 0)

get_diff_plot =  function(feature1, feature2, neg_control, NAS_data, xlab, ylab, swap = FALSE, big_only = FALSE, return_p = TRUE, return_plot = FALSE, title=NULL, title_left=NULL, hline = NULL, ylims = NULL, breaks = NULL, labels = NULL, scale_y = NULL) {
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
  
  n_t <- get_n_t_value_line(NAS_data, feature1, feature2)
  colnames(n_t) <- c("true_diff")
  n_t$ID <- seq.int(nrow(n_t))
  
  plot <- ggplot()
  
  if(!is.null(hline)) {
    print(hline)
    plot <- plot +
      geom_hline(yintercept = hline, lty=2, col="#aaaaaa") +
      geom_hline(yintercept = -hline, lty=2, col="#aaaaaa") +
      geom_rect(aes(xmin = 0, xmax=max(n_t$ID), ymin=-hline, ymax=hline), col="grey", alpha=0.1, linetype=0)
  }
  
  plot <- plot + 
    geom_line(aes(x = n_t$ID, y = n_t$true_diff), col="RoyalBlue") + 
    labs(x = xlab, y=ylab) + 
    theme(
      # axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      panel.background = element_rect( fill = "#f5f5f5"),
      panel.grid.major = element_line(colour="#ffffff"),
      panel.grid.minor = element_line(colour="#ffffff")
    )
  
  plot <- plot + scale_x_continuous(expand = c(0,0), labels = c())
  
  if(!is.null(ylims) | !is.null(labels) | !is.null(breaks)) {
    plot <- plot + scale_y_continuous(limits = ylims, labels = labels, breaks = breaks)
  }
  if(!is.null(scale_y)) {
    plot <- plot +
      coord_cartesian(ylim = scale_y)
  }
  plot <- title_styling(plot, title, title_left)
  return(plot)
}




get_neg_control_plot <- function(NAS_data, neg_control, feature1, feature2, limit_y = NULL, order=FALSE) {
  library(ggplot2)
  library(reshape2)
  true_diffs <- get_true_diffs(NAS_data, "PSI_no_PTC", "PSI_het_PTC")
  sims_diffs <- get_sims_diffs(NAS_data, "PSI_no_PTC", "PSI_het_PTC", true_diffs)
  
  sims_diffs <- melt(sims_diffs)
  
  sims_diffs = sims_diffs[!is.na(sims_diffs$value),]
  sims_diffs$id <- numextract(sims_diffs$variable)
  
  # print(head(sims_diffs))
  
  uniques <- length(unique(sims_diffs$id))
  # print(length(uniques))
  l <- list()
  for(i in 1:uniques) {
    l[[i]] = c()
  }
  
  # sims <- data.frame(matrix(NA, nrow = max_diff_length, ncol = length(neg_control)))
  for(i in 1:nrow(sims_diffs)) {
    row <- sims_diffs[i,]
    id = row$id
    l[[id]] <- c(l[[id]], row$value)
    # print(id)
    
    # do stuff with row
  }
  positives <- list()
  for (i in 1:length(l)) {
    entry = l[[i]]
    pos <- sum(entry >= 0)
    positives[[i]] = pos
  }
  
  sims_diffs <- mutate(sims_diffs, positives = NA, order = NA)
  sims_diffs$id <- as.numeric(sims_diffs$id)
  
  # print(head(sims_diffs))
  
  for (i in 1:nrow(sims_diffs)) {
    row <- sims_diffs[i,]
    # print(as.numeric(row$id))
    # print(positives[as.numeric(row$id)])
    # print(positives[as.numeric(row$id)])
    sims_diffs[i, "positives"] = positives[row$id]
  }
  
  sims_diffs <- arrange(sims_diffs, desc(positives))
  rank = 0
  current_id = 0
  for (i in 1:nrow(sims_diffs)) {
    row <- sims_diffs[i,]
    if (row$id != current_id) {
      current_id = row$id
      rank = rank + 1
    }
    sims_diffs[i, "order"] = rank
  }
  
  if(order) {
    xcat = "order"
  } else {
    xcat = "id"
  }
  
  plot <- ggplot() +
    geom_point(aes(x = sims_diffs[[xcat]], y = sims_diffs$value), size=0.2, colour=ifelse(sims_diffs$value > 0, "RoyalBlue", "red")) +
    geom_hline(yintercept = 0, lty=2) +
    labs(x="PTC", y = "True PSIdiff - simulant PSIdiff") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "#f5f5f5"),
      panel.grid.major = element_line(colour="#f5f5f5"),
      panel.grid.minor = element_line(colour="#ffffff")
    )
  
  if(!is.null(limit_y)) {
    plot <- plot + coord_cartesian(ylim=limit_y)
  }
  
  return(plot)
}


get_n_t =  function(feature1, feature2, neg_control, NAS_data, title, swap = FALSE, big_only = FALSE, return_p = TRUE, return_plot = FALSE) {
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
  n_t = rep(NA, length(neg_control))
  for (index in 1:length(neg_control)) {
    true_diff = NAS_data[index, feature1] - NAS_data[index, feature2]
    if (abs(true_diff) >= threshold) {
      sim_diffs = neg_control[[index]][feature1] - neg_control[[index]][feature2]
      sim_diffs = sim_diffs[!is.na(sim_diffs)]
      n = sum((true_diff > sim_diffs) & (true_diff != sim_diffs))
      t = length(sim_diffs[true_diff != sim_diffs])
      n_t[index] = n/t      
    }
  }  
  plot <- hist(n_t, col = "RoyalBlue", breaks = 20, main = title, xlab = "Proportion of simulants that show an effect smaller or equal to true effect.")
  plot
  abline(v = 0.5, lty = 2, lwd = 2)
  print(sum((n_t > 0.5) & (!is.na(n_t))))
  print(sum(!is.na(n_t)))
  p = binom.test(sum((n_t > 0.5) & (!is.na(n_t))), sum(!is.na(n_t)), alternative = "g")$p.value
  text(0.7, y_pos, label = paste("p ~ ", round(p, 3), sep = ""), cex = 2)
  if (return_p == TRUE) {
    return(p)
  } else if (return_plot == TRUE) {
    return(plot)
  } else {
    return(n_t)
  }
}

get_n_t_plot =  function(feature1, feature2, neg_control, NAS_data, xlab, swap = FALSE, big_only = FALSE, return_p = TRUE, return_plot = FALSE, title=NULL, title_left=NULL) {
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
  
  library(ggplot2)
  plot <- ggplot() + 
    geom_histogram(aes(n_t),binwidth = 0.02,fill = chosen_colour,color = "black") + 
    geom_vline(xintercept = 0.5, lty=2) +
    coord_cartesian(xlim = c(0,1)) +
    labs(x = xlab, y="Frequency") + 
    plot_theme(font_size = 12)
  
  plot <- title_styling(plot, title, title_left)
  return(plot)
}

get_n_t_value <- function(neg_control, feature1, feature2, threshold) {
  n_t = rep(NA, length(neg_control))
  for (index in 1:length(neg_control)) {
    # get the true difference
    true_diff = NAS_data[index, feature1] - NAS_data[index, feature2]
    # if the true difference is bigger than the threshold
    if (abs(true_diff) >= threshold) {
      # get the differences for each of the simulations
      sim_diffs = neg_control[[index]][feature1] - neg_control[[index]][feature2]
      sim_diffs = sim_diffs[!is.na(sim_diffs)]
      n = sum((true_diff > sim_diffs) & (true_diff != sim_diffs))
      t = length(sim_diffs[true_diff != sim_diffs])
      n_t[index] = n/t      
    }
  }
  print(sum((n_t > 0.5) & (!is.na(n_t))))
  print(sum(!is.na(n_t)))
  p = binom.test(sum((n_t > 0.5) & (!is.na(n_t))), sum(!is.na(n_t)), alternative = "g")$p.value
  print(p)
  # return(n_t)
}

get_n_t_value_line <- function(NAS_data, feature1, feature2) {
  n_t = rep(NA, length(NAS_data))
  diffs <- data.frame(true = double())
  for (index in 1:length(neg_control)) {
    # get the true difference
    true_diff = NAS_data[index, feature1] - NAS_data[index, feature2]
    diffs <- rbind(diffs, data.frame(true_diff))
  }
  return(diffs)
}


get_n_t_visual =  function(feature1, feature2, neg_control, NAS_data, swap = FALSE) {
  par(mfrow = c(4, 4))
  if (swap == TRUE) {
    temp = feature1
    feature1 = feature2
    feature2 = temp
  }
  for (index in 1:length(neg_control)) {
    true_diff = NAS_data[index, feature1] - NAS_data[index, feature2]
    sim_diffs = neg_control[[index]][feature1] - neg_control[[index]][feature2]
    sim_diffs = sim_diffs[!is.na(sim_diffs)]
    x_min = min(sim_diffs, true_diff)
    x_max = max(sim_diffs, true_diff)
    hist(sim_diffs, xlim = c(x_min, x_max), breaks = 20, col = "RoyalBlue", main = "")
    abline(v = true_diff, lwd = 2, col = "orange")
  }
}

get_sims_diffs <- function(NAS_data, feature1, feature2, true_diffs) {
  max_diff_length = 0
  for (index in 1:length(neg_control)) {
    sim_diffs = neg_control[[index]][feature1] - neg_control[[index]][feature2]
    sim_diffs = sim_diffs[!is.na(sim_diffs)]
    if(length(sim_diffs) > max_diff_length) {
      max_diff_length = length(sim_diffs)
    }
  }
  
  sims <- data.frame(matrix(NA, nrow = max_diff_length, ncol = length(neg_control)))
  for (index in 1:length(neg_control)) {
    sim_diffs = neg_control[[index]][feature1] - neg_control[[index]][feature2]
    sim_diffs = sim_diffs[!is.na(sim_diffs)]
    true_diff = true_diffs[index, "true_diff"]
    sim_diffs = true_diff - sim_diffs
    length(sim_diffs) <- max_diff_length
    col_id = paste("X", index, sep="")
    sims[col_id] <- sim_diffs
  }
  return(sims)
}


get_true_diffs <- function(NAS_data, feature1, feature2) {
  diffs <- data.frame(true = double())
  for (index in 1:length(neg_control)) {
    # get the true difference
    true_diff = NAS_data[index, feature1] - NAS_data[index, feature2]
    diffs <- rbind(diffs, data.frame(true_diff))
  }
  return(diffs)
}


individual_changes_plot <- function(NAS_data, feature1, feature2, threshold, xlab, ylab, col1  = NULL, col2 = NULL, xlab1=NULL, xlab2=NULL) {
  dots <- data.frame(xpos = double(), ypos = double(), diff = double())
  lines <- data.frame(x0 = double(), x1 = double(), y0 = double(), y1 = double())
  for (exon in 1:dim(NAS_data)[1]) {
    exon_name = NAS_data[exon, 1]
    no_PTC = NAS_data[exon, feature1]
    het_PTC = NAS_data[exon, feature2]
    x_pos = 0
    x_pos2 = 1
    if (is.finite(no_PTC) && is.finite(het_PTC)) {
      if (abs(no_PTC - het_PTC) > threshold) {
        diff = no_PTC - het_PTC
        dots[nrow(dots) + 1,] = list(x_pos,no_PTC,diff)
        dots[nrow(dots) + 1,] = list(x_pos2,het_PTC,diff)
        lines[nrow(lines) + 1, ] = list(x_pos, x_pos2, no_PTC, het_PTC)
      }
    }
  }
  
  if(is.null(col1)) {
    col1 = "RoyalBlue"
  }
  if(is.null(col2)) {
    col2 = "red"
  }
  
  if(is.null(xlab1) || is.null(xlab2)) {
    labels <- c("PTC -/-", "PTC-/+")
  } else {
    labels <- c(xlab1, xlab2)
  }
  
  plot <- ggplot() + 
    geom_point(aes(x = dots$xpos, y = dots$ypos), size=0.8, colour = ifelse(dots$diff > 0, col1, col2))
  plot <- plot +
    geom_segment(data = lines, aes(x = x0, y = y0, xend = x1, yend = y1), colour = ifelse(lines$y1 > lines$y0, col2, col1)) +
    scale_x_continuous(breaks=c(0, 1), labels=labels) +
    labs(x = xlab, y = ylab) +
    theme(
      panel.background = element_rect( fill = "#f5f5f5"),
      panel.grid.major = element_line(colour="#ffffff"),
      panel.grid.minor = element_line(colour="#ffffff")
    )
  
  return(plot)
}

individual_changes_pvals_plot <- function(NAS_data) {
  pvals <- c()
  for (i in seq(0, 10, 0.1)) {
    changes_PSI_binom_test = plot_individual_change(NAS_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "", "PSI", i, 100, big_changes_binom_test = TRUE)
    pvals <- c(pvals, changes_PSI_binom_test$p.value)
  }
  
  max(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC)
  x <- seq(0, 10, 0.1)
  y <- log(pvals)
  
  plot <- ggplot() +
    geom_point(aes(x= x, y =y), col = ifelse(y < log(0.05), "RoyalBlue", "red")) + 
    geom_hline(yintercept = log(0.05), lty=2) +
    scale_x_continuous(name = "PSIdiff threshold", breaks=seq(0, 10, by = 1), labels = seq(0, 10, by = 1)) +
    labs(y="Bionmial test log P value") + 
    theme(
      # axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      panel.background = element_rect( fill = "#f5f5f5"),
      panel.grid.major = element_line(colour="#ffffff"),
      panel.grid.minor = element_line(colour="#ffffff")
    )
  return(plot)
}


numextract <- function(string){ 
  library(stringr)
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 


overlap_significance = function(overlap, pool_size, n1, n2, sim_number) {
  sim_overlaps = vector()
  pool = seq(1:pool_size)
  for (sim in 1:sim_number) {
    sample1 = sample(pool, n1)
    sample2 = sample(pool, n2)
    sim_overlaps = c(sim_overlaps, length(intersect(sample1, sample2)))
  }
  greater = length(sim_overlaps[sim_overlaps > overlap])
  p = (greater + 1)/(sim_number + 1)
  print(greater)
  print(p)
}

perform_tests = function(NAS_data) {
  print("H: Heterozygotes have lower PSI than PTC- homozygotes")
  print(wilcox.test(NAS_data$PSI_het_PTC, NAS_data$PSI_no_PTC, alternative = "l", paired = TRUE)$p.value)
  print("H: PTC+ homozygotes have lower PSI than heterozygotes")
  print(wilcox.test(NAS_data$PSI_w_PTC, NAS_data$PSI_het_PTC, alternative = "l", paired = TRUE)$p.value)
  print("H: PTC+ homozygotes have lower PSI than PTC- homozygotes")
  print(wilcox.test(NAS_data$PSI_w_PTC, NAS_data$PSI_no_PTC, alternative = "l", paired = TRUE)$p.value)
  #TBA: a repeated measures ANOVA type test
  
  print("H: Heterozygotes have higher normalized read count than PTC- homozygotes")
  print(wilcox.test(NAS_data$norm_count_het_PTC, NAS_data$norm_count_no_PTC, alternative = "g", paired = TRUE)$p.value)
  print("H: PTC+ homozygotes have higher normalized read count than heterozygotes")
  print(wilcox.test(NAS_data$norm_count_w_PTC, NAS_data$norm_count_het_PTC, alternative = "g", paired = TRUE)$p.value)
  print("H: PTC+ homozygotes have higher normalized read count than PTC- homozygotes")
  print(wilcox.test(NAS_data$norm_count_w_PTC, NAS_data$norm_count_no_PTC, alternative = "g", paired = TRUE)$p.value)
}

plot_diff_hist_het <- function(title=NULL, title_left=NULL) {
  library(ggplot2)
  plot <- ggplot() + 
    geom_histogram(aes(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC),binwidth = 1,fill = chosen_colour,color = "black") + 
    scale_x_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
    labs(x = "PSI PTC-/- - PTC-/+", y="Frequency") + 
    plot_theme(font_size = 12)
  
  plot <- title_styling(plot, title, title_left)
  return(plot)
}


plot_diff_hist_het_rpm <- function(feature1, feature2, title=NULL, title_left=NULL, binwidth=NULL, xlab=NULL, ylab=NULL) {
  if(is.null(binwidth)) {
    bindwidth <- 1
  }
  if(is.null(xlab)) {
    xlab <- ""
  }
  if(is.null(ylab)) {
    ylab <- ""
  }
  
  library(ggplot2)
  plot <- ggplot() + 
    geom_histogram(aes(NAS_data[[feature1]] - NAS_data[[feature2]]), binwidth = binwidth,fill = chosen_colour,color = "black") + 
    # scale_x_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
    coord_cartesian(ylim=c(0,20)) +
    labs(x = xlab, y= ylab) + 
    plot_theme(font_size = 12) +
    geom_vline(xintercept = 0.025, col="red") +
    geom_vline(xintercept = -0.025, col="red")
  
  plot <- title_styling(plot, title, title_left)
  return(plot)
}

plot_diff_hist_het_zoom <- function(title=NULL, title_left=NULL) {
  library(ggplot2)
  plot <- ggplot() + 
    geom_histogram(aes(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC),binwidth = 1,fill = chosen_colour,color = "black") + 
    coord_cartesian(ylim=c(0, 50)) +
    scale_x_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
    labs(x = "PSI PTC-/- - PTC-/+", y="Frequency") + 
    plot_theme(font_size = 12)
  
  plot <- title_styling(plot, title, title_left)
  return(plot)
}

plot_diff_hists_het = function(NAS_data, title) {
  if (!is.null(NAS_data$PSI_het_PTC) & !is.null(NAS_data$PSI_no_PTC)) {
    hist(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, breaks = 100, col = chosen_colour, main = title, xlab = "PSI(-/-) - PSI(-/+)")
  }
}

plot_diff_hists_homo = function(NAS_data, title) {
  if (!is.null(NAS_data$PSI_het_PTC) & !is.null(NAS_data$PSI_w_PTC)) {
    hist(NAS_data$PSI_het_PTC - NAS_data$PSI_w_PTC, breaks = 100, col = chosen_colour, main = title, xlab = "PSI(-/+) - PSI(+/+)")
  }
}

plot_diff_RPM_hists_het = function(NAS_data, title) {
  if (!is.null(NAS_data$norm_count_het_PTC) & !is.null(NAS_data$norm_count_no_PTC)) {
    hist(NAS_data$norm_count_no_PTC - NAS_data$norm_count_het_PTC, breaks = 100, col = chosen_colour, main = title, xlab = "RPMskip(-/-) - RPMskip(-/+)")
  }
}

plot_diff_RPM_hists_homo = function(NAS_data, title) {
  if (!is.null(NAS_data$norm_count_het_PTC) & !is.null(NAS_data$norm_count_w_PTC)) {
    hist(NAS_data$norm_count_het_PTC - NAS_data$norm_count_w_PTC, breaks = 100, col = chosen_colour, main = title, xlab = "RPMskip(-/+) - RPMskip(+/+)")
  }
}

plot_individual_change = function(NAS_data, w_PTC_name, het_PTC_name, no_PTC_name, title, ylab, threshold, max_value, reverse = FALSE, big_changes_binom_test = FALSE) {
  big_changes = vector()
  colours = c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  plot.new()
  plot.window(c(0,2), c(0, max_value))
  axis(1, at = c(0, 1, 2), labels = c("PTC-/-", "PTC-/+", "PTC+/+"))
  axis(2, at = seq(0, max_value, by = max_value/10), labels = seq(0, max_value, by = max_value/10))
  title(main = title, xlab = "genotype", ylab = ylab) 
  counter = 1
  coin_toss_total = 0
  coin_toss_greater = 0
  for (exon in 1:dim(NAS_data)[1]) {
    change_colour = FALSE
    exon_name = NAS_data[exon, 1]
    no_PTC = NAS_data[exon, no_PTC_name]
    het_PTC = NAS_data[exon, het_PTC_name]
    w_PTC = NAS_data[exon, w_PTC_name]
    x_pos = jitter(0)
    x_pos2 = jitter(1)
    x_pos3 = jitter(2)
    colour = colours[counter]
    if (is.finite(no_PTC) && is.finite(het_PTC)) {
      if (abs(no_PTC - het_PTC) > threshold) {
        points(x_pos, no_PTC, pch = 19, col = colour)
        points(x_pos2, het_PTC, pch = 19, col = colour)
        segments(x0 = x_pos, x1 = x_pos2, y0 = no_PTC, y1 = het_PTC, col = colour, lwd = 1.5)
        change_colour = TRUE
        if (is.finite(w_PTC)) {
          points(x_pos3, w_PTC, pch = 19, col = colour)
          segments(x0 = x_pos2, x1 = x_pos3, y0 = het_PTC, y1 = w_PTC, col = colour, lwd = 1.5)
        }
        coin_toss_total = coin_toss_total + 1
        if (no_PTC > het_PTC) {
          if (!reverse) {
            big_changes = c(big_changes, exon_name)
          }
          coin_toss_greater = coin_toss_greater + 1
        }
        else {
          if (reverse) {
            big_changes = c(big_changes, exon_name)
          }
        }
      }
    }
      if (is.finite(w_PTC) && is.finite(het_PTC)) {
        if (abs(w_PTC - het_PTC) > threshold) {
          points(x_pos2, het_PTC, pch = 19, col = colour)
          points(x_pos3, w_PTC, pch = 19, col = colour)
          segments(x0 = x_pos2, x1 = x_pos3, y0 = het_PTC, y1 = w_PTC, col = colour, lwd = 1.5)
          change_colour = TRUE
          if (is.finite(no_PTC)) {
            points(x_pos, no_PTC, pch = 19, col = colour)
            segments(x0 = x_pos, x1 = x_pos2, y0 = no_PTC, y1 = het_PTC, col = colour, lwd = 1.5)
          }
        }
      }
    if (change_colour == TRUE) {
      counter = counter + 1
      if (counter > length(colours)) {
        counter = 1
      }
    }
  }
  print(coin_toss_greater)
  print(coin_toss_total)
  if (reverse == TRUE) {
    alternative = "less"
  }
  else {
    alternative = "greater"
  }
  
  btest <- binom.test(x = coin_toss_greater, n = coin_toss_total, p = 0.5, alternative = alternative)

  if(big_changes_binom_test) {
    return(btest)
  } else {
    return(big_changes) 
  }
}

plot_individual_change_lim = function(NAS_data, w_PTC_name, het_PTC_name, no_PTC_name, title, ylab, threshold, max_value, reverse = FALSE, big_changes_binom_test = FALSE) {
  big_changes = vector()
  colours = c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  plot.new()
  plot.window(c(0,1), c(0, max_value))
  axis(1, at = c(0, 1), labels = c("PTC-/-", "PTC-/+"))
  axis(2, at = seq(0, max_value, by = max_value/10), labels = seq(0, max_value, by = max_value/10))
  title(main = title, xlab = "genotype", ylab = ylab) 
  counter = 1
  coin_toss_total = 0
  coin_toss_greater = 0
  for (exon in 1:dim(NAS_data)[1]) {
    change_colour = FALSE
    exon_name = NAS_data[exon, 1]
    no_PTC = NAS_data[exon, no_PTC_name]
    het_PTC = NAS_data[exon, het_PTC_name]
    x_pos = jitter(0)
    x_pos2 = jitter(1)
    colour = colours[counter]
    if (is.finite(no_PTC) && is.finite(het_PTC)) {
      if (abs(no_PTC - het_PTC) > threshold) {
        diff = no_PTC - het_PTC
        points(x_pos, no_PTC, pch = 19, col = ifelse(diff > 0, "blue", "red"))
        points(x_pos2, het_PTC, pch = 19, col = ifelse(diff > 0, "blue", "red"))
        segments(x0 = x_pos, x1 = x_pos2, y0 = no_PTC, y1 = het_PTC, col = ifelse(diff > 0, "blue", "red"), lwd = 1.5)
        change_colour = TRUE
        # if (is.finite(w_PTC)) {
        #   points(x_pos3, w_PTC, pch = 19, col = colour)
        #   segments(x0 = x_pos2, x1 = x_pos3, y0 = het_PTC, y1 = w_PTC, col = colour, lwd = 1.5)
        # }
        coin_toss_total = coin_toss_total + 1
        if (no_PTC > het_PTC) {
          if (!reverse) {
            big_changes = c(big_changes, exon_name)
          }
          coin_toss_greater = coin_toss_greater + 1
        }
        else {
          if (reverse) {
            big_changes = c(big_changes, exon_name)
          }
        }
      }
    }
    # if (is.finite(w_PTC) && is.finite(het_PTC)) {
    #   if (abs(w_PTC - het_PTC) > threshold) {
    #     points(x_pos2, het_PTC, pch = 19, col = colour)
    #     points(x_pos3, w_PTC, pch = 19, col = colour)
    #     segments(x0 = x_pos2, x1 = x_pos3, y0 = het_PTC, y1 = w_PTC, col = colour, lwd = 1.5)
    #     change_colour = TRUE
    #     if (is.finite(no_PTC)) {
    #       points(x_pos, no_PTC, pch = 19, col = colour)
    #       segments(x0 = x_pos, x1 = x_pos2, y0 = no_PTC, y1 = het_PTC, col = colour, lwd = 1.5)
    #     }
    #   }
    # }
    if (change_colour == TRUE) {
      counter = counter + 1
      if (counter > length(colours)) {
        counter = 1
      }
    }
  }
  print(coin_toss_greater)
  print(coin_toss_total)
  if (reverse == TRUE) {
    alternative = "less"
  }
  else {
    alternative = "greater"
  }
  
  btest <- binom.test(x = coin_toss_greater, n = coin_toss_total, p = 0.5, alternative = alternative)
  
  if(big_changes_binom_test) {
    return(btest)
  } else {
    return(big_changes) 
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

# plot_individual_change_lim(NAS_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "", "", 5, 100, reverse = T, big_changes_binom_test = F)

prepare_controls <- function(set) {
  set <- transform(
    set,
    PSI_no_PTC = PSI_no_PTC*100,
    PSI_w_PTC = PSI_w_PTC*100,
    PSI_het_PTC = PSI_het_PTC*100,
    norm_count_w_PTC = norm_count_w_PTC * 1000000,
    norm_count_het_PTC = norm_count_het_PTC * 1000000,
    norm_count_no_PTC = norm_count_no_PTC * 1000000,
    norm_count_w_PTC_incl = norm_count_w_PTC_incl * 1000000,
    norm_count_het_PTC_incl = norm_count_het_PTC_incl * 1000000,
    norm_count_no_PTC_incl = norm_count_no_PTC_incl * 1000000
  )
}

prepare_dataset = function(datafile) {
  NAS_data = read.csv(datafile, stringsAsFactors = FALSE, sep = "\t")
  max_sample_count = max(NAS_data$sample_count)
  NAS_data = NAS_data[NAS_data$ptc_count > 0 & NAS_data$sample_count > (0.5 * max_sample_count), ]
  NAS_data$PSI_w_PTC = NAS_data$PSI_w_PTC * 100
  NAS_data$PSI_het_PTC = NAS_data$PSI_het_PTC * 100
  NAS_data$PSI_no_PTC = NAS_data$PSI_no_PTC * 100
  NAS_data$norm_count_w_PTC = NAS_data$norm_count_w_PTC * 1000000
  NAS_data$norm_count_het_PTC = NAS_data$norm_count_het_PTC * 1000000
  NAS_data$norm_count_no_PTC = NAS_data$norm_count_no_PTC * 1000000
  if ("norm_count_w_PTC_incl" %in% colnames(NAS_data)) {
    NAS_data$norm_count_w_PTC_incl = NAS_data$norm_count_w_PTC_incl * 1000000
    NAS_data$norm_count_het_PTC_incl = NAS_data$norm_count_het_PTC_incl * 1000000
    NAS_data$norm_count_no_PTC_incl = NAS_data$norm_count_no_PTC_incl * 1000000
  }
  NAS_data = NAS_data[order(NAS_data$id),]
  return(NAS_data)
}

prepare_neg_control <- function(out_list) {
  neg_control <- lapply(out_list, prepare_controls)
  return(neg_control)
}

read_in_simulations = function(simulant_prefix, simulant_number, ids, column_names) {
  sim_mat = matrix(NA, 1, 13)
  colnames(sim_mat) = column_names
  out_list = rep(list(sim_mat), length(ids))
  for (sim in 1:simulant_number) {
    print(sim)
    current_data = read.csv(paste(simulant_prefix, sim, ".txt", sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    for (exon in 1:dim(current_data)[1]) {
      current_id = current_data[exon, "id"]
      if (current_id %in% ids) {
        out_list[[match(current_id, ids)]] = rbind(out_list[[match(current_id, ids)]], current_data[exon, ])
      }
    }
  }
  
  prepared_out_list <- prepare_neg_control(out_list)
  return(prepared_out_list)
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




####### tests

chosen_colour = "RoyalBlue"


#PSI
NAS_data = prepare_dataset("results/clean_run_2/clean_run__analysis_final_output.txt")
NAS_data_shift = prepare_dataset("results/clean_run_2/clean_run_out_of_frame__analysis_final_output.txt")
neg_control = read_in_simulations("results/clean_run_2/simulation_output/final_output_simulation_", 100, NAS_data$id, colnames(NAS_data))

# perform tests between genotypes
perform_tests(NAS_data)

# samples with a PTC+/+ sample
nrow(NAS_data[NAS_data$PSI_w_PTC != 'NaN',])

# summary of data
summary(NAS_data)

# medians of data
median(NAS_data$PSI_het_PTC)
median(NAS_data$PSI_no_PTC)

# ask whether there is a significant difference between PSI for non-PTC homozygotes and PTC heterozygotes
wilcox.test(NAS_data$PSI_het_PTC, NAS_data$PSI_no_PTC, paired=T)
wilcox.test(NAS_data$PSI_het_PTC, NAS_data$PSI_no_PTC, paired=T, alternative = "less")
wilcox.test(NAS_data$PSI_het_PTC, NAS_data$PSI_no_PTC, paired=T, alternative = "greater")

# ask whether there is a significant difference between RPMinclude for het PTC and homo no PTC
wilcox.test(NAS_data$norm_count_het_PTC_incl, NAS_data$norm_count_no_PTC, paired=T, alternative = "greater")
pdf('results/graphs/rpminclude_het_ptc_homo_no_ptc.pdf', width=10)
boxplot(NAS_data$norm_count_het_PTC_incl, NAS_data$norm_count_no_PTC_incl, ylim=c(0, 5), labels = c("PTC -/+", "PTC-/-"))
dev.off()

# ask whether there is a significant difference between RPMskip for her PTC and homo no PTC
wilcox.test(NAS_data$norm_count_het_PTC, NAS_data$norm_count_no_PTC, paired=T, alternative="greater")

#histogram of psi ptc-/- - ptc-/+
hist_psi_no_ptc_minus_psi_het_ptc <- plot_diff_hist_het(title="A", title_left=T)
ggsave("results/graphs/hist_psi_diff_no_ptc_het_ptc.pdf", plot=hist_psi_no_ptc_minus_psi_het_ptc, width=8, height=5)

#histogram of psi ptc-/- - ptc-/+ with zoom
psi_diff_plot <- get_diff_plot("PSI_no_PTC", "PSI_het_PTC", neg_control, NAS_data, xlab="PTC", ylab="PSIdiff", hline = 5, ylims = c(-20, 100), breaks = seq(-20, 100, 10), labels = seq(-20, 100, 10))
ggsave("results/graphs/psi_diff.pdf", plot = psi_diff_plot, width=8, height = 5 )
# hist_psi_no_ptc_minus_psi_het_ptc_title <- plot_diff_hist_het(title="A", title_left=T)
# hist_psi_no_ptc_minus_psi_het_ptc_zoom_title <- plot_diff_hist_het_zoom(title="B", title_left=T)
# hist_psi_no_het_with_zoom <- grid.arrange(hist_psi_no_ptc_minus_psi_het_ptc_title, hist_psi_no_ptc_minus_psi_het_ptc_zoom_title, ncol=2)
# ggsave("results/graphs/hist_psi_no_het_with_zoom.pdf", plot=hist_psi_no_het_with_zoom, width=12, height=5)

# histogram of the simulants with psi difference less than real ptc-/- - ptc-/+
simulants_less_real_ptc_no_minus_het <- get_n_t_plot("PSI_no_PTC", "PSI_het_PTC", neg_control, NAS_data, xlab="Proportion of simulants with PSI for pseudoPTC-/- - pseudoPTC-/+\nless than or equal to PTC-/- - PTC-/+", swap = FALSE, title="A", title_left=T)
# same for rpmskip
simulants_less_real_ptc_no_minus_het_rpm_skip <- get_n_t_plot("norm_count_no_PTC", "norm_count_het_PTC", neg_control, NAS_data, xlab="Proportion of simulants with RPMskip for pseudoPTC-/- - pseudoPTC-/+\nless than or equal to PTC-/- - PTC-/+",  swap = T, title="B", title_left =T)


p <- get_n_t_plot("PSI_no_PTC", "PSI_het_PTC", neg_control, NAS_data, xlab="Proportion of simulants with PSI for pseudoPTC-/- - pseudoPTC-/+\nless than or equal to PTC-/- - PTC-/+", swap = FALSE, title="A", title_left=T, return_plot = F, return_p)
p

get_neg_control_plot(NAS_data, neg_control, "PSI_no_PTC", "PSI_het_PTC")

# get the number of ptcs where the real PTC diff if greater than 50% of simulants
get_n_t_value(neg_control, "PSI_no_PTC", "PSI_het_PTC", 0)


# just psi
ggsave('results/graphs/simulants_less_real_ptc_no_minus_het.pdf', plot = simulants_less_real_ptc_no_minus_het, width=8, height=5)
# both
sim_plot <- grid.arrange(simulants_less_real_ptc_no_minus_het, simulants_less_real_ptc_no_minus_het_rpm_skip, ncol=2)
ggsave('results/graphs/simulantions_less_than_real_diff.pdf', plot = sim_plot, width=14, height=5)

# big changes in psi

individuals_large_effects <- individual_changes_plot(NAS_data, "PSI_no_PTC", "PSI_het_PTC", 5, "Genotype", "PSI")
individuals_large_effects

big_changes_psi = plot_individual_change(NAS_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "", "", 5, 100, reverse = F, big_changes_binom_test = TRUE)
big_changes_psi
big_changes_psi = plot_individual_change(NAS_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "", "", 5, 100, reverse = T, big_changes_binom_test = F)
big_changes_RPMskip = plot_individual_change(NAS_data, "norm_count_w_PTC", "norm_count_het_PTC", "norm_count_no_PTC", "Exons with >0.025 change between any two categories", "RPMskip", 0.025, 5, reverse = TRUE)


# combine plots
library(ggpubr)
plot <- ggarrange(psi_diff_plot, individuals_large_effects, widths = c(3, 2), ncol = 2, nrow = 1, labels = c("A", "B"))
ggsave('results/graphs/psi_diff_large_effect_plot.pdf', plot = plot, width=12, height=7)


# is rpmskip greater for PTC-/- variants?
wilcox.test(NAS_data$norm_count_no_PTC, NAS_data$norm_count_het_PTC, paired = T, alternative = "greater")

# rpmskip compared with simulants
get_n_t("norm_count_no_PTC", "norm_count_het_PTC", neg_control, NAS_data, "title", swap = T, big_only = FALSE, return_p = T, return_plot = FALSE)

# diff plots
rpmskip_plot <- get_diff_plot("norm_count_no_PTC", "norm_count_het_PTC", neg_control, NAS_data, xlab="PTC", ylab="RPMskipDiff", hline = 0.025, ylims = c(-1.5, 1), breaks = seq(-1.5, 1, 0.25), labels = seq(-1.5, 1, 0.25), scale_y = c(-0.5, 0.5))
rpmskip_big_changes <- individual_changes_plot(NAS_data, "norm_count_no_PTC", "norm_count_het_PTC", 0.025, "Genotype", "RPMskip", col="red", col2="RoyalBlue")
rpmincl_plot <- get_diff_plot("norm_count_no_PTC_incl", "norm_count_het_PTC_incl", neg_control, NAS_data, xlab="PTC", ylab="RPMinclDiff", hline = 0, ylims = c(-25, 40), breaks = seq(-25, 40, 5), labels = seq(-25, 40, 5), scale_y = c(-10, 10))
plot <- ggarrange(
  ggarrange(rpmincl_plot, rpmskip_plot, nrow=2, labels=c("", "B")),
  ncol=2,
  rpmskip_big_changes,
  labels=c("A", "C")
)
ggsave('results/graphs/rpm_diff_large_effect_plot.pdf', plot = plot, width=12, height=10)

# get the number of PTCs for which the real RPMskip is greater than 50% of simulants
get_n_t("norm_count_no_PTC", "norm_count_het_PTC", neg_control, NAS_data, "title", swap=T)


# cases for PTCs shifted by 1 nucleotide
# get the number in the correction direction
big_changes_psi_shift = plot_individual_change(NAS_data_shift, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "", "", 5, 100, reverse = F, big_changes_binom_test = TRUE)
big_changes_RPMskip_shift = plot_individual_change(NAS_data_shift, "norm_count_w_PTC", "norm_count_het_PTC", "norm_count_no_PTC", "Exons with >0.025 change between any two categories", "RPMskip", 0.0375, 5, reverse = TRUE, big_changes_binom_test = TRUE)
big_changes_psi_shift
big_changes_RPMskip_shift

# show individual large cases
shiftptc_psi_changes <- individual_changes_plot(NAS_data_shift, "PSI_no_PTC", "PSI_het_PTC", 5, "Genotype", "PSI", xlab1 = "shiftPTC-/-", xlab2="shiftPTC-/+")
shiftptc_rpmskip_changes <- individual_changes_plot(NAS_data_shift, "norm_count_no_PTC", "norm_count_het_PTC", 0.0375, "Genotype", "RPMskip", xlab1 = "shiftPTC-/-", xlab2="shiftPTC-/+", col1="red", col2="RoyalBlue")
plot <- ggarrange(shiftptc_psi_changes, shiftptc_rpmskip_changes, ncol=2, labels=c("A", "B"))
ggsave("results/graphs/shiftptc_individual_changes.pdf", plot = plot, width=10, height=7)



neg_control_plot <- get_neg_control_plot(NAS_data, neg_control, "PSI_no_PTC", "PSI_het_PTC")
neg_control_plot_zoom <- get_neg_control_plot(NAS_data, neg_control, "PSI_no_PTC", "PSI_het_PTC", limit = c(-3, 3))
plot <- ggarrange(neg_control_plot, neg_control_plot_zoom, widths = c(2, 2), ncol = 2, nrow = 1, labels = c("A", "B"))
neg_control_plot_zoom
ggsave('results/graphs/simulants_psi_diff.pdf', plot = plot, width=12, height=7)
ggsave('results/graphs/simulants_psi_diff.eps', plot = plot, width=12, height=7)


# plot showing the p values for each percentage cutoff
plot = individual_changes_pvals_plot(NAS_data)
ggsave("results/graphs/psi_difference_binomial_test_p_values.eps", plot = plot, width=10, height=7)  



#analysis of the exons that show a big change
#is the overlap between the ones that show a big change in RPMskip and those that show a big change
#in PSI greater than you'd expect by chance?
overlap = intersect(big_changes_PSI, big_changes_RPMskip)
overlap
overlap_significance(length(overlap), dim(NAS_data)[1], length(big_changes_PSI), length(big_changes_RPMskip), 10000)
print(sort(overlap))

expression = read.csv("results/clean_run_2/clean_run_FANTOM_expression_per_transcript.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
expression_analysis(overlap, NAS_data$exon, expression)


# RS
#are the ones that show a big change expressed at lower levels than those that don't (Geuvadis data)?
all_exons_expression = NAS_data$norm_count_no_PTC_incl[!(NAS_data$exon %in% overlap)] + NAS_data$norm_count_no_PTC[!(NAS_data$exon %in% overlap)]
overlap_expression = NAS_data$norm_count_no_PTC_incl[NAS_data$exon %in% overlap] + NAS_data$norm_count_no_PTC[NAS_data$exon %in% overlap]
wilcox.test(all_exons_expression, overlap_expression, alt = "g")

#are the ones that show a big change less affected by NMD?
NMD_effect = NAS_data$norm_count_no_PTC_incl - NAS_data$norm_count_het_PTC_incl
all_exons_NMD = NMD_effect[!(NAS_data$exon %in% overlap)]
overlap_NMD = NMD_effect[NAS_data$exon %in% overlap]
wilcox.test(all_exons_NMD, overlap_NMD, alt = "g")










