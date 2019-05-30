library(reshape2)
library(ggplot2)

p_from_z <- function(z) {
  return(2*pnorm(-abs(z)))
}

get_region_z <- function(data, set) {
  real = data[data$simulation == "real",]
  sims = data[data$simulation != "real",]
  z0_2 <- (real$X0.2 - mean(sims$X0.2)) / sd(sims$X0.2)
  z3_69 <- (real$X3.69 - mean(sims$X3.69)) / sd(sims$X3.69)
  z70 <- (real$X70. - mean(sims$X70.)) / sd(sims$X70.)
  
  zs <- data.frame(set, z0_2, z3_69, z70)
  colnames(zs) <- c("set", "z_0_2","z_3_69", "z_70")
  return(zs)
}

get_region_hit_z <- function(data) {
  real = data[data$simulation == "real",]
  sims = data[data$simulation != "real",]
  z <- (real$ese_hit_count - mean(sims$ese_hit_count)) / sd(sims$ese_hit_count)
  zs <- data.frame(z)
  colnames(zs) <- c("region hits")
  return(zs)
}


get_region_ps <- function(data, set) {
  real = data[data$simulation == "real",]
  sims = data[data$simulation != "real",]
  p0_2 <- (nrow(sims[sims$X0.2 >= real$X0.2,]) + 1)/(nrow(sims) + 1)
  p3_69 <- (nrow(sims[sims$X3.69 >= real$X3.69,]) + 1)/(nrow(sims) + 1)
  p70 <- (nrow(sims[sims$X70. >= real$X70.,]) + 1)/(nrow(sims) + 1)
  # z0_2 <- (real$X0.2 - mean(sims$X0.2)) / sd(sims$X0.2)
  # z3_69 <- (real$X3.69 - mean(sims$X3.69)) / sd(sims$X3.69)
  # z70 <- (real$X70. - mean(sims$X70.)) / sd(sims$X70.)
  
  ps <- data.frame(set, p0_2, p3_69, p70)
  colnames(ps) <- c("set", "p_0_2","p_3_69", "p_70")
  return(ps)
}

get_region_hit_p <- function(data) {
  real = data[data$simulation == "real",]
  sims = data[data$simulation != "real",]
  # z <- (real$ese_hit_count - mean(sims$ese_hit_count)) / sd(sims$ese_hit_count)
  # p <- data.frame(p_from_z(z))
  
  head(real)
  p = (nrow(sims[sims$ese_hit_count >= real$ese_hit_count,]) + 1) / (nrow(sims) + 1)
  print(p)
  # colnames(p) <- c("region p")
  # return(p)p = (nrow(sims[sims$ese_hit_count >= real$ese_hit_count,]) + 1) / (nrow(sims) + 1)
  # colnames(p) <- c("region p")
  return(p)
}

###

locations_file <- "results/clinvar2/clinvar_ptc_location_simulation.csv"
file <- read.csv(locations_file, head=T)
clinvar_regions_zs <- get_region_z(file, "ClinVar")
clinvar_regions_ps <- get_region_ps(file, "ClinVar")
clinvar_regions_zs
clinvar_regions_ps

real = file[file$simulation == "real",]
sims = file[file$simulation != "real",]


int3_file <- "results/clinvar2/clinvar_ptc_location_simulation_CaceresHurstESEs_INT3_ese_overlaps.csv"
rescue_file <- "results/clinvar2/clinvar_ptc_location_simulation_RESCUE_eses_ese_overlaps.csv"
ke400_file <- "results/clinvar2/clinvar_ptc_location_simulation_Ke400_eses_ese_overlaps.csv"

int3 <- read.csv(int3_file, head=T)
rescue <- read.csv(rescue_file, head=T)
ke400 <- read.csv(ke400_file, head=T)

z_int3 <- get_region_z(int3, "INT3")
z_rescue <- get_region_z(rescue, "RESCUE")
z_ke400 <- get_region_z(ke400, "Ke400")


get_region_ps(int3, "INT3")
get_region_ps(rescue, "RESCUE")
get_region_ps(ke400, "Ke400")


ese_hits_int3_file <- "results/clinvar2/clinvar_ese_hit_simulation_3_69_CaceresHurstESEs_INT3.csv"
ese_hits_rescue_file <- "results/clinvar2/clinvar_ese_hit_simulation_3_69_RESCUE_eses.csv"
ese_hits_ke400_file <- "results/clinvar2/clinvar_ese_hit_simulation_3_69_Ke400_eses.csv"
ese_hits_int3 <- read.csv(ese_hits_int3_file, head=T)
ese_hits_rescue <- read.csv(ese_hits_rescue_file, head=T)
ese_hits_ke400 <- read.csv(ese_hits_ke400_file, head=T)

get_region_hit_z(ese_hits_int3)
get_region_hit_z(ese_hits_rescue)
get_region_hit_z(ese_hits_ke400)
get_region_hit_p(ese_hits_int3)
get_region_hit_p(ese_hits_rescue)
get_region_hit_p(ese_hits_ke400)


ese_hits <- rbind(z_int3, z_rescue, z_ke400)
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


kgenomes_location_file <- read.csv('results/clinvar2/1000_genomes_simulations.csv', head=T)
kgenomes_regions_zs <- get_region_z(kgenomes_location_file, "1000 Genomes")

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

ggsave('results/graphs/nonsense_mutations_locations_simulations.eps', plot=plot, width=10, height=7)



