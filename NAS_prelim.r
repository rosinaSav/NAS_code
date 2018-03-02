library(vioplot)

my_vioplot = function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                       horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                       lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                       at, add = FALSE, wex = 1, drawRect = TRUE, cexAxis = 1, mtext = NULL, ylab = "") 
{ datas <- list(x, ...)
  n <- length(datas)
  if (length(col) == 1) {
    col <- rep(col, n)}
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25, na.rm = TRUE)
    q3[i] <- quantile(data, 0.75, na.rm = TRUE)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label, cex.axis = cexAxis)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  if (!missing(mtext)) {
    text_pos1 = at[seq(1, length(at), by = 2)]
    text_pos2 = at[seq(2, length(at), by = 2)]
    text_pos = (text_pos1 + text_pos2)/2
    mtext(mtext, side = 3, line = 0, at = text_pos)
  }
  mtext(ylab, side = 2, line = 2.5)
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}

prepare_dataset = function(datafile) {
  NAS_data = read.csv(datafile, stringsAsFactors = FALSE, sep = "\t")
  max_sample_count = max(NAS_data$sample_count)
  NAS_data = NAS_data[NAS_data$ptc_count > 0 & NAS_data$sample_count > (0.5 * max_sample_count), ]
  NAS_data$norm_count_w_PTC = NAS_data$norm_count_w_PTC * 1000000
  NAS_data$norm_count_het_PTC = NAS_data$norm_count_het_PTC * 1000000
  NAS_data$norm_count_no_PTC = NAS_data$norm_count_no_PTC * 1000000
  return(NAS_data)
}

NAS_data = prepare_dataset("results/clean_run/clean_run_final_output.txt")
max_sample_count = max(NAS_data$sample_count)
NAS_data = NAS_data[NAS_data$ptc_count > 0 & NAS_data$sample_count > (0.5 * max_sample_count), ]
NAS_data$norm_count_w_PTC = NAS_data$norm_count_w_PTC * 1000000
NAS_data$norm_count_het_PTC = NAS_data$norm_count_het_PTC * 1000000
NAS_data$norm_count_no_PTC = NAS_data$norm_count_no_PTC * 1000000
wilcox.test(NAS_data$PSI_het_PTC, NAS_data$PSI_no_PTC, alternative = "t", paired = TRUE)
wilcox.test(NAS_data$norm_count_het_PTC, NAS_data$norm_count_no_PTC, alternative = "t", paired = TRUE)
wilcox.test(NAS_data$PSI_w_PTC, NAS_data$PSI_het_PTC, alternative = "t", paired = TRUE)
wilcox.test(NAS_data$norm_count_w_PTC, NAS_data$norm_count_het_PTC, alternative = "t", paired = TRUE)

summary(NAS_data)

par(mfrow = c(2, 1))
boxplot(NAS_data$PSI_w_PTC, NAS_data$PSI_het_PTC, NAS_data$PSI_no_PTC, col = "limegreen", main = "PSI", names = c("hmz PTC", "htz", "hmz no PTC"))
boxplot(NAS_data$norm_count_w_PTC, NAS_data$norm_count_het_PTC, NAS_data$norm_count_no_PTC, col = "RoyalBlue", main = "normalized read count", names = c("hmz PTC", "htz", "hmz no PTC"))

my_vioplot(NAS_data$PSI_w_PTC, NAS_data$PSI_het_PTC, NAS_data$PSI_no_PTC, names = c("hmz PTC", "htz", "hmz no PTC"), col = "limegreen", ylab = "PSI") 
my_vioplot(NAS_data$norm_count_w_PTC, NAS_data$norm_count_het_PTC, NAS_data$norm_count_no_PTC, names = c("hmz PTC", "htz", "hmz no PTC"), col = "limegreen", ylab = "normalized read count") 

NAS_ESEs_data = read.csv("results/clean_run/clean_run_CaceresHurstESEs_INT3_final_output.txt", stringsAsFactors = FALSE, sep = "\t")
NAS_ESEs_data = NAS_ESEs_data[NAS_ESEs_data$ptc_count > 0 & NAS_ESEs_data$sample_count > (0.5 * max_sample_count), ]
