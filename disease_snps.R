p_from_z <- function(z) {
  return(2*pnorm(-abs(z)))
}

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

get_region_hit_z <- function(file) {
  print(head(file))
  real = file[file$simulation == "real",]
  sims = file[file$simulation != "real",]
  z <- (real$both_ese_hit_count - mean(sims$both_ese_hit_count)) / sd(sims$both_ese_hit_count)
  zs <- data.frame(z)
  colnames(zs) <- c("region hits")
  return(zs)
}

file <- read.csv('results/clinvar/clinvar_ese_hit_simulation_3_69.csv', head=T)
get_region_hit_z(file)


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

locations_file <- "results/clinvar/clinvar_simulations.csv"
clinvar_regions_zs <- get_region_z(locations_file, "ClinVar")
clinvar_regions_ps <- get_region_ps(locations_file, "ClinVar")
clinvar_regions_zs
clinvar_regions_ps

int3 <- "results/clinvar/clinvar_simulations_ese_overlaps.csv"
int3_clinvar_ese_hit_zs <- get_region_z(int3, "INT3")
int3_clinvar_ese_hit_ps <- get_region_ps(int3, "INT3")
int3_clinvar_ese_hit_zs
int3_clinvar_ese_hit_ps

rescue <- "results/clinvar/clinvar_simulations_RESCUE_eses_overlaps.csv"
rescue_clinvar_ese_hit_zs <- get_region_z(rescue, "RESCUE")
rescue_clinvar_ese_hit_ps <- get_region_ps(rescue, "RESCUE")
rescue_clinvar_ese_hit_zs
rescue_clinvar_ese_hit_ps

ke <- "results/clinvar/clinvar_simulations_Ke400_eses_overlaps.csv"
ke_clinvar_ese_hit_zs <- get_region_z(ke, "Ke400")
ke_clinvar_ese_hit_ps <- get_region_ps(ke, "Ke400")
ke_clinvar_ese_hit_zs
ke_clinvar_ese_hit_ps

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
plot
ggsave('results/graphs/nonsense_mutations_ese_hits_locations_simulations.eps', plot=plot, width=10, height=7)



file <- read.csv('results/clinvar/clinvar_ese_hit_simulation_3_69.csv', head=T)
get_region_hit_z(file)

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



