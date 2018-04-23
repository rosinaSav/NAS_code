compare_the_two_measures = function(NAS_data, title) {
  cor_test = cor.test(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, NAS_data$norm_count_no_PTC - NAS_data$norm_count_het_PTC, method = "spearman")
  plot(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, NAS_data$norm_count_no_PTC - NAS_data$norm_count_het_PTC, pch = 19, col = "RoyalBlue", xlab = "PSI(-/-) - PSI(-/+)", ylab = "RPMskip(-/-) - RPMskip(-/+)", main = title)
  abline(h = 0)
  abline(v = 0)
  text(x = 10, y = 1, labels = paste("rho ~ ", round(cor_test$estimate, 3), "\np ~ ", round(cor_test$p.value, 3), sep = ""), adj = 0, cex = 1.3)
}

get_n_t =  function(feature1, feature2, neg_control, NAS_data, title, swap = FALSE, big_only = FALSE) {
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
  hist(n_t, col = "RoyalBlue", breaks = 20, main = title, xlab = "Proportion of simulants that show an effect smaller or equal to true effect.")
  abline(v = 0.5, lty = 2, lwd = 2)
  p = binom.test(sum((n_t > 0.5) & (!is.na(n_t))), sum(!is.na(n_t)), alternative = "g")$p.value
  text(0.7, y_pos, label = paste("p ~ ", round(p, 3), sep = ""), cex = 2)
  return(p)
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
  if ("norm_count_w_PTC_incl" %in% colnames(NAS_data)) {
    NAS_data$norm_count_w_PTC_incl = NAS_data$norm_count_w_PTC_incl * 1000000
    NAS_data$norm_count_het_PTC_incl = NAS_data$norm_count_het_PTC_incl * 1000000
    NAS_data$norm_count_no_PTC_incl = NAS_data$norm_count_no_PTC_incl * 1000000
  }
  NAS_data = NAS_data[order(NAS_data$id),]
  return(NAS_data)
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
  return(out_list)
}

ese_overlap_simulation_plot <- function(snp_file, monomorphic_file) {
  real <- snp_file[snp_file$id == "real",]
  sims_snp <- snp_file[snp_file$id != "real",]
  sims_mono <- monomorphic_file[monomorphic_file$id != "real",]
  
  real$snps_overlapping_eses <- real$snps_overlapping_eses/1179
  sims_snp$snps_overlapping_eses <- sims_snp$snps_overlapping_eses/1179
  sims_mono$snps_overlapping_eses <- sims_mono$snps_overlapping_eses/1179

  max_overlaps <- max(real$snps_overlapping_eses, sims_mono$snps_overlapping_eses, sims_snp$snps_overlapping_eses)
  min_overlaps <- min(real$snps_overlapping_eses, sims_mono$snps_overlapping_eses, sims_snp$snps_overlapping_eses)

  hist(sims_mono$snps_overlapping_eses, freq=FALSE, breaks=seq(0, 0.25, 0.0025), xlim=c(0,0.21), ylim=c(0,120), col="RoyalBlue", xlab="Proportion of cases overlapping > 1 ESE motif", main="")
  hist(sims_snp$snps_overlapping_eses, breaks=seq(0, 1, 0.0025), add=T, col="orange", freq=FALSE)
  abline(v = real$snps_overlapping_eses, lty=2, lwd=2)
  legend(0.04,120,
         legend=c(paste("Allele frequency/ancestral allele matched SNPs (n=", nrow(sims_snp), ")", sep=""), paste("Ancestral allele matched monomorphic sites (n=", nrow(sims_mono), ")", sep=""), "Real PTC causing SNPs"),
         col=c("orange", "RoyalBlue", "black"),
         lty=c(NA,NA,2),
         lwd=c(NA,NA,2),
         cex=0.8,
         box.lty=0,
         pch=c(15,15,NA)
  )
}
ese_overlap_simulation_plot(snp_sim_file, mono_sim_file)

ese_overlap_simulation_tests <- function(snp_file, monomorphic_file) {
  real <- snp_file[snp_file$id == "real",]
  sims_snp <- snp_file[snp_file$id != "real",]
  sims_mono <- monomorphic_file[monomorphic_file$id != "real",]

  p_snp_sim <- nrow(sims_snp[sims_snp$snps_overlapping_eses > real$snps_overlapping_eses,]) / nrow(sims_snp)
  p_mono_sim <- nrow(sims_mono[sims_mono$snps_overlapping_eses < real$snps_overlapping_eses,]) / nrow(sims_mono)

  print("H: PTCs overlap ESEs significantly more than randomly selected non-PTC SNPs.")
  print(p_snp_sim)
  print("H: PTCs overlap ESEs significantly less than randomly selected non-SNP sites.")
  print(p_mono_sim)
}



chosen_colour = "RoyalBlue"

#PSI
NAS_data = prepare_dataset("results/clean_run/clean_run_final_output.txt")
neg_control = read_in_simulations("results/clean_run/simulation_output/final_output_simulation_", 30, NAS_data$id, colnames(NAS_data))
perform_tests(NAS_data)
summary(NAS_data)

par(mfrow = c(2, 2))
p_PSI = get_n_t("PSI_no_PTC", "PSI_het_PTC", neg_control, NAS_data, "PSI-/- - PSI-/+", swap = FALSE)
p_RPMskip = get_n_t("norm_count_no_PTC", "norm_count_het_PTC", neg_control, NAS_data, "RPMskip-/- - RPMskip-/+", swap = TRUE)
p_PSI_big = get_n_t("PSI_no_PTC", "PSI_het_PTC", neg_control, NAS_data, "PSI-/- - PSI-/+ (5+)", swap = FALSE, big_only = 5)
p_RPMskip_big = get_n_t("norm_count_no_PTC", "norm_count_het_PTC", neg_control, NAS_data, "RPMskip-/- - RPMskip-/+ (0.05+)", swap = TRUE, big_only = 0.05)

get_n_t_visual("PSI_no_PTC", "PSI_het_PTC", neg_control, NAS_data, swap = FALSE)
get_n_t_visual("norm_count_no_PTC", "norm_count_het_PTC", neg_control, NAS_data, swap = TRUE)

par(mfrow = c(1, 3))
hist(NAS_data$PSI_no_PTC, main = "PTC-/-", xlab = "PSI", col = chosen_colour, breaks = 20)
hist(NAS_data$PSI_het_PTC, main = "PTC-/+", xlab = "PSI", col = chosen_colour, breaks = 20)
hist(NAS_data$PSI_w_PTC, main = "PTC+/+", xlab = "PSI", col = chosen_colour, breaks = 20)

jpeg("results/clean_run/minusminus_vs_het.jpeg", width = 23.4, height = 15.4, units = 'cm', res = 300)
par(mfrow = c(5, 6))
par("mar" = c(3.1, 2.1, 2.1, 1))
plot_diff_hists_het(NAS_data, "PTC-/PTC- - PTC-/PTC+")
for (sim in 1:30) {
  plot_diff_hists_het(neg_control[[sim]], sim)
}
dev.off()

jpeg("results/clean_run/het_vs_plusplus.jpeg", width = 24.4, height = 15.4, units = 'cm', res = 300)
par(mfrow = c(5, 6))
par("mar" = c(3.1, 2.1, 2.1, 1))
plot_diff_hists_homo(NAS_data, "PTC-/PTC+ - PTC+/PTC+")
for (sim in 1:30) {
  plot_diff_hists_homo(neg_control[[sim]], sim)
}
dev.off()

par(mfrow = c(1,2))
hist(NAS_data$PSI_no_PTC - NAS_data$PSI_het_PTC, breaks = 100, col = chosen_colour, main = "Is PSI higher for PTC-/PTC- than for PTC-/PTC+?", xlab = "PSI(-/-) - PSI(-/+)")
hist(NAS_data$PSI_het_PTC - NAS_data$PSI_w_PTC, breaks = 100, col = chosen_colour, main = "Is PSI higher for PTC-/PTC+ than for PTC+/PTC+?", xlab = "PSI(-/+) - PSI(+/+)")

graphics.off()
plot_individual_change(NAS_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "Exons with >5% change between any two categories", "PSI", 5, 100)

#RPMskip
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

jpeg("results/clean_run/minusminus_vs_het_RPM.jpeg", width = 25.4, height = 15.4, units = 'cm', res = 300)
par(mfrow = c(5, 6))
par("mar" = c(3.1, 2.1, 2.1, 1))
plot_diff_RPM_hists_het(NAS_data, "-/-hom - het (RPM)")
for (sim in 1:30) {
  plot_diff_RPM_hists_het(neg_control[[sim]], sim)
}
dev.off()

jpeg("results/clean_run/het_vs_plusplus_RPM.jpeg", width = 25.4, height = 15.4, units = 'cm', res = 300)
par(mfrow = c(5, 6))
par("mar" = c(3.1, 2.1, 2.1, 1))
plot_diff_RPM_hists_homo(NAS_data, "het - +/+hom (RPM)")
for (sim in 1:30) {
  plot_diff_hists_homo(neg_control[[sim]], sim)
}
dev.off()
plot_individual_change(NAS_data, "norm_count_w_PTC", "norm_count_het_PTC", "norm_count_no_PTC", "Exons with >0.05 change between any two categories", "RPMskip", 0.05, 6, reverse = TRUE)
plot_individual_change(neg_control, "norm_count_w_PTC", "norm_count_het_PTC", "norm_count_no_PTC", "Exons with >0.02 change between any two categories (negative control)", "RPMskip", 0.02, 1.12, reverse = TRUE)

#ESE analysis
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

par(mfrow = c(2, 2))
plot_individual_change(NAS_ESEs_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "in ESEs", "PSI", 5, 100)
plot_individual_change(NAS_no_ESEs_data, "PSI_w_PTC", "PSI_het_PTC", "PSI_no_PTC", "not in ESEs", "PSI", 5, 100)
plot_individual_change(NAS_ESEs_data, "norm_count_w_PTC", "norm_count_het_PTC", "norm_count_no_PTC", "in ESEs", "RPMskip", 0.05, 2, reverse = TRUE)
plot_individual_change(NAS_no_ESEs_data, "norm_count_w_PTC", "norm_count_het_PTC", "norm_count_no_PTC", "not in ESEs", "RPMskip", 0.05, 7, reverse = TRUE)

#RPMinclude
par(mfrow = c(1,1))
hist(NAS_data$norm_count_no_PTC_incl - NAS_data$norm_count_het_PTC_incl, breaks = 100, col = chosen_colour, main = "Is RPM(incl) higher for PTC-/PTC- than for PTC-/PTC+?", xlab = "RPMincl(-/-) - RPMincl(-/+)")
wilcox.test(NAS_data$norm_count_no_PTC_incl, NAS_data$norm_count_het_PTC_incl, alternative = "greater", paired = TRUE)

# ESE overlap simulations
snp_sim_file <- read.csv('results/clean_run_2/ese_overlap_simulation/snp_simulation/ese_overlap_snp_simulation.csv', head=T)
mono_sim_file <- read.csv('results/clean_run_2/ese_overlap_simulation/monomorphic_sim/ese_overlap_monomorphic_simulation.csv', head=T)
jpeg('results/clean_run_2/plots/ese_overlap_simulation.jpg', width=25, height = 15, units="cm", res=300)
ese_overlap_simulation_plot(snp_sim_file, mono_sim_file)
dev.off()

ese_overlap_simulation_tests(snp_sim_file, mono_sim_file)
