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
perform_tests(NAS_data)
summary(NAS_data)

NAS_ESEs_data = prepare_dataset("results/clean_run/clean_run_CaceresHurstESEs_INT3_final_output.txt")
perform_tests(NAS_ESEs_data)

NAS_no_ESEs_data = prepare_dataset("results/clean_run/clean_run_CaceresHurstESEs_INT3_complement_final_output.txt")
perform_tests(NAS_no_ESEs_data)
