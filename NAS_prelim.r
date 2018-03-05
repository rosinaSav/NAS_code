compare_the_two_measures = function(NAS_data, title) {
  cor_test = cor.test(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, NAS_data$norm_count_no_PTC - NAS_data$norm_count_het_PTC, method = "spearman")
  plot(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, NAS_data$norm_count_no_PTC - NAS_data$norm_count_het_PTC, pch = 19, col = "RoyalBlue", xlab = "PSI(-/-) - PSI(-/+)", ylab = "RPMskip(-/-) - RPMskip(-/+)", main = title)
  abline(h = 0)
  abline(v = 0)
  text(x = 10, y = 1, labels = paste("rho ~ ", round(cor_test$estimate, 3), "\np ~ ", round(cor_test$p.value, 3), sep = ""), adj = 0, cex = 1.3)
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
  print(wilcox.test(NAS_data$norm_count_het_PTC, NAS_data$norm_count_no_PTC, alternative = "l", paired = TRUE)$p.value)
  print("H: PTC+ homozygotes have higher normalized read count than heterozygotes")
  print(wilcox.test(NAS_data$norm_count_w_PTC, NAS_data$norm_count_het_PTC, alternative = "l", paired = TRUE)$p.value)
  print("H: PTC+ homozygotes have higher normalized read count than PTC- homozygotes")
  print(wilcox.test(NAS_data$norm_count_w_PTC, NAS_data$norm_count_no_PTC, alternative = "l", paired = TRUE)$p.value)
}

plot_individual_change = function(NAS_data, w_PTC_name, het_PTC_name, no_PTC_name, title, ylab, threshold, max_value, reverse = FALSE) {
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
          coin_toss_greater = coin_toss_greater + 1
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
  print(binom.test(x = coin_toss_greater, n = coin_toss_total, p = 0.5, alternative = alternative))
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
  return(NAS_data)
}

chosen_colour = "RoyalBlue"

#PSI
NAS_data = prepare_dataset("results/clean_run/clean_run_final_output.txt")
neg_control = prepare_dataset("results/clean_run/clean_run_simulate_ptc_snps_bam_analysis/final_output_simulation_1.txt")
perform_tests(NAS_data)
perform_tests(neg_control)
summary(NAS_data)
summary(neg_control)
par(mfrow = c(2, 3))
hist(NAS_data$PSI_no_PTC, main = "PTC-/-", xlab = "PSI", col = chosen_colour, breaks = 20)
hist(NAS_data$PSI_het_PTC, main = "PTC-/+", xlab = "PSI", col = chosen_colour, breaks = 20)
hist(NAS_data$PSI_w_PTC, main = "PTC+/+", xlab = "PSI", col = chosen_colour, breaks = 20)
hist(neg_control$PSI_no_PTC, main = "PTC-/- (neg. control)", xlab = "PSI", col = chosen_colour, breaks = 20)
hist(neg_control$PSI_het_PTC, main = "PTC-/+ (neg. control)", xlab = "PSI", col = chosen_colour, breaks = 20)
hist(neg_control$PSI_w_PTC, main = "PTC+/+ (neg.control)", xlab = "PSI", col = chosen_colour, breaks = 20)

par(mfrow = c(2,2))
hist(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, breaks = 100, col = chosen_colour, main = "Is PSI higher for PTC-/PTC- than for PTC-/PTC+?", xlab = "PSI(-/-) - PSI(-/+)")
hist(NAS_data$PSI_het_PTC - NAS_data$PSI_w_PTC, breaks = 100, col = chosen_colour, main = "Is PSI higher for PTC-/PTC+ than for PTC+/PTC+?", xlab = "PSI(-/+) - PSI(+/+)")
hist(neg_control$PSI_no_PTC - neg_control$PSI_het_PTC, breaks = 100, col = chosen_colour, main = "Negative control", xlab = "PSI(-/-) - PSI(-/+)")
hist(neg_control$PSI_het_PTC - neg_control$PSI_w_PTC, breaks = 100, col = chosen_colour, main = "Negative control", xlab = "PSI(-/+) - PSI(+/+)")

plot_individual_change(NAS_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "Exons with >5% change between any two categories", "PSI", 5, 100)
plot_individual_change(neg_control, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "Negative control", "PSI", 5, 100)

#normalized read counts
par(mfrow = c(1, 2))
compare_the_two_measures(NAS_data, "true data")
compare_the_two_measures(neg_control, "neg. control")

par(mfrow = c(2, 3))
hist(NAS_data$norm_count_no_PTC, main = "PTC-/-", xlab = "RPM that support skipping", col = chosen_colour, breaks = 20)
hist(NAS_data$norm_count_het_PTC, main = "PTC-/+", xlab = "RPM that support skipping", col = chosen_colour, breaks = 20)
hist(NAS_data$norm_count_w_PTC, main = "PTC+/+", xlab = "RPM that support skipping", col = chosen_colour, breaks = 20)
hist(neg_control$norm_count_no_PTC, main = "PTC-/- (neg. control)", xlab = "RPM that support skipping", col = chosen_colour, breaks = 20)
hist(neg_control$norm_count_het_PTC, main = "PTC-/+ (neg. control)", xlab = "RPM that support skipping", col = chosen_colour, breaks = 20)
hist(neg_control$norm_count_w_PTC, main = "PTC+/+ (neg.control)", xlab = "RPM that support skipping", col = chosen_colour, breaks = 20)

par(mfrow = c(2,2))
hist(NAS_data$norm_count_no_PTC - NAS_data$norm_count_het_PTC, breaks = 100, col = chosen_colour, main = "Is RPM(skip) lower for PTC-/PTC- than for PTC-/PTC+?", xlab = "RPMskip(-/-) - RPMskip(-/+)")
hist(NAS_data$norm_count_het_PTC - NAS_data$norm_count_w_PTC, breaks = 100, col = chosen_colour, main = "Is RPM(skip) lower for PTC-/PTC+ than for PTC+/PTC+?", xlab = "RPMskip(-/+) - RPMskip(+/+)")
hist(neg_control$norm_count_no_PTC - neg_control$norm_count_het_PTC, breaks = 100, col = chosen_colour, main = "Negative control", xlab = "RPMskip(-/-) - RPMskip(-/+)")
hist(neg_control$norm_count_het_PTC - neg_control$norm_count_w_PTC, breaks = 100, col = chosen_colour, main = "Negative control", xlab = "RPMskip(-/+) - RPMskip(+/+)")

plot_individual_change(NAS_data, "norm_count_w_PTC", "norm_count_het_PTC", "norm_count_no_PTC", "Exons with >0.1 change between any two categories", "RPMskip", 0.05, 6, reverse = TRUE)
plot_individual_change(neg_control, "norm_count_w_PTC", "norm_count_het_PTC", "norm_count_no_PTC", "Negative control", "RPMskip", 0.05, 8, reverse = TRUE)

NAS_ESEs_data = prepare_dataset("results/clean_run/clean_run_CaceresHurstESEs_INT3_final_output.txt")
perform_tests(NAS_ESEs_data)
NAS_no_ESEs_data = prepare_dataset("results/clean_run/clean_run_CaceresHurstESEs_INT3_complement_final_output.txt")
perform_tests(NAS_no_ESEs_data)
wilcox.test(NAS_ESEs_data$PSI_no_PTC - NAS_ESEs_data$PSI_het_PTC, NAS_no_ESEs_data$PSI_no_PTC - NAS_no_ESEs_data$PSI_het_PTC, alternative = "greater")
wilcox.test(NAS_ESEs_data$norm_count_no_PTC - NAS_ESEs_data$norm_count_het_PTC, NAS_no_ESEs_data$norm_count_no_PTC - NAS_no_ESEs_data$norm_count_het_PTC, alternative = "less")

par(mfrow = c(2, 2))
hist(NAS_ESEs_data$PSI_no_PTC - NAS_ESEs_data$PSI_het_PTC, col = chosen_colour, breaks = 20, main = "PSI, in ESEs", xlab = "PSI(-/-) - PSI(-/+)")
hist(NAS_no_ESEs_data$PSI_no_PTC - NAS_no_ESEs_data$PSI_het_PTC, col = chosen_colour, breaks = 20, main = "PSI, not in ESEs", xlab = "PSI(-/-) - PSI(-/+)")
hist(NAS_ESEs_data$norm_count_no_PTC - NAS_ESEs_data$norm_count_het_PTC, col = chosen_colour, breaks = 20, main = "RPMskip, in ESEs", xlab = "RPMskip(-/-) - RPMskip(-/+)")
hist(NAS_no_ESEs_data$norm_count_no_PTC - NAS_no_ESEs_data$norm_count_het_PTC, col = chosen_colour, breaks = 20, main = "RPMskip, not in ESEs", xlab = "RPMskip(-/-) - RPMskip(-/+)")

par(mfrow = c(1, 2))
plot_individual_change(NAS_ESEs_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "Exons with >5% change between any two categories (in ESEs)", "PSI", 5, 100)
plot_individual_change(NAS_no_ESEs_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "Exons with >5% change between any two categories (not in ESEs)", "PSI", 5, 100)
