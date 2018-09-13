p_from_z <- function(z) {
  return(2*pnorm(-abs(z)))
}

info_file <- read.table('results/clinvar/disease__analysis_ptc_info.txt', head=T, sep=",")
nrow(info_file[info_file$min_dist >= 3 & info_file$min_dist <= 69,])
nrow(info_file[info_file$min_dist >= 3 & info_file$min_dist <= 69,]) / nrow(info_file)

get_region_z <- function(filepath, set) {
  file <- read.csv(filepath, head=T)
  real = file[file$sim_id == "real",]
  sims = file[file$sim_id != "real",]
  z0_2 <- (real$all_0.2.bp - mean(sims$all_0.2.bp)) / sd(sims$all_0.2.bp)
  z3_69 <- (real$all_3.69.bp - mean(sims$all_3.69.bp)) / sd(sims$all_3.69.bp)
  z70 <- (real$all_70..bp - mean(sims$all_70..bp)) / sd(sims$all_70..bp)
  
  zs <- data.frame(set, z0_2, z3_69, z70)
  colnames(zs) <- c("set", "z_0_2","z_3_69", "z_70")
  return(zs)
}

get_region_ps <- function(filepath, set) {
  file <- read.csv(filepath, head=T)
  real = file[file$sim_id == "real",]
  sims = file[file$sim_id != "real",]
  z0_2 <- (real$all_0.2.bp - mean(sims$all_0.2.bp)) / sd(sims$all_0.2.bp)
  z3_69 <- (real$all_3.69.bp - mean(sims$all_3.69.bp)) / sd(sims$all_3.69.bp)
  z70 <- (real$all_70..bp - mean(sims$all_70..bp)) / sd(sims$all_70..bp)
  
  ps <- data.frame(set, p_from_z(z0_2), p_from_z(z3_69), p_from_z(z70))
  colnames(ps) <- c("set", "z_0_2","z_3_69", "z_70")
  return(ps)
}

clinvar_regions_zs <- get_region_z("results/clinvar/clinvar_simulations.csv", "ClinVar")
clinvar_regions_ps <- get_region_ps("results/clinvar/clinvar_simulations.csv", "ClinVar")
clinvar_regions_ps

kgenomes_regions_zs <- get_region_z("results/clinvar/1000_genomes_simulations.csv", "1000 Genomes")



nonsense_mutation_regions <- rbind(kgenomes_regions_zs, clinvar_regions_zs)
nonsense_mutation_regions_melt <- melt(nonsense_mutation_regions)

plot <- ggplot(nonsense_mutation_regions_melt, aes(x=variable, fill=set, y=value)) +   
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Exon region (nucleotides)", y="Z") +
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +
  geom_hline(yintercept = 1.96, lty=2) +
  geom_hline(yintercept = -1.96, lty=2) +
  geom_hline(yintercept = 0, lty=1) +
  scale_x_discrete(labels = c("0-2", "3-69", "70+")) + 
  theme(legend.title=element_blank())

ggsave('results/graphs/nonsense_mutations_locations_simulations.pdf', plot=plot, width=10, height=7)


int3_clinvar_ese_hit_zs <- get_region_z("results/clinvar/clinvar_simulations_ese_overlaps.csv", "INT3")
rescue_clinvar_ese_hit_zs <- get_region_z("results/clinvar/clinvar_simulations_RESCUE_eses_overlaps.csv", "RESCUE")
ke_clinvar_ese_hit_zs <- get_region_z("results/clinvar/clinvar_simulations_Ke400_eses_overlaps.csv", "Ke400")

ese_hits <- rbind(int3_clinvar_ese_hit_zs, rescue_clinvar_ese_hit_zs, ke_clinvar_ese_hit_zs)
ese_hits_melt <- melt(ese_hits)

plot <- ggplot(ese_hits_melt, aes(x=variable, fill=set, y=value)) +   
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Exon region (nucleotides)", y="Z") +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) +
  geom_hline(yintercept = 1.96, lty=2) +
  geom_hline(yintercept = -1.96, lty=2) +
  geom_hline(yintercept = 0, lty=1) +
  scale_x_discrete(labels = c("0-2", "3-69", "70+")) +
  theme(legend.title=element_blank())
ggsave('results/graphs/nonsense_mutations_ese_hits_locations_simulations.pdf', plot=plot, width=10, height=7)

