file1 <- read.table("./results/disease_snps/disease__analysis_ptc_info.txt", head=T, sep=",")
hist(file1$min_dist, breaks=seq(min(file1$min_dist),max(file1$min_dist)+10, by=10), freq=F)

file2 <- read.table("./results/disease_snps/disease__analysis_other_info.txt", head=T, sep=",")
hist(file2$min_dist, breaks=seq(min(file2$min_dist),max(file2$min_dist)+10, by=10), freq = F)


nrow(file2[file2$min_dist<70,])/nrow(file2)


real_prop <- nrow(file1[file1$min_dist<70,])/nrow(file1)

sim_props <- c()
for (i in 1:100) {
  print(i)
  file <- read.table(paste('results/disease_snps/ptc_location_simulation/disease_location_simulation_', i, '.bed', sep=''), sep="\t")
  prop <- nrow(file[file$V19 < 50,])/nrow(file)
  sim_props <- c(sim_props, prop)
}
sim_props

greater <- sum(sim_props > real_prop)
p <- (greater + 1) / (length(sim_props) + 1)
p


file <- read.table('results/clinvar/ptc_distances_disease.txt', head=T, sep="\t")
head(file)

disease = file[file$disease == 1,]
others = file[file$disease == 0,]

nrow(others[others$exon_dist <= 3,]) / nrow(others)
nrow(disease[disease$exon_dist <= 3,]) / nrow(disease)

nrow(others[others$exon_dist > 3 & others$exon_dist <= 69,]) / nrow(others)
nrow(disease[disease$exon_dist > 3 & disease$exon_dist <= 69,]) / nrow(disease)


wilcox.test(file$exon_dist[file$disease ==1], file$exon_dist[file$disease==0])

nrow(file[file$disease ==1,])

# are disease ptcs in genes with more exons?
# absolute number of exons in disease genes is higher than non disease
wilcox.test(disease$total_exons, others$total_exons)

# but disease genes tend to be longer, does this make a difference
# if we are comparing between ptcs?
wilcox.test(disease$cds_length, others$cds_length)

median(disease$total_exons)
median(others$total_exons)


# does the position in the gene have an influence on where ptcs occur?
first_half = nrow(file[file$gene_left_length < file$gene_right_length,])
binom.test(x = first_half, n = nrow(file), p = 0.5, alternative = "g")

first_half_disease = nrow(file[file$gene_left_length < file$gene_right_length & file$disease == 1 ,])
binom.test(x = first_half_disease, n = nrow(file[file$disease == 1,]), p = 0.5, alternative = "g")

first_half_non_disease = nrow(file[file$gene_left_length < file$gene_right_length & file$disease != 1 ,])
binom.test(x = first_half_non_disease, n = nrow(file[file$disease != 1,]), p = 0.5, alternative = "g")


file <- read.table("results/clinvar/compare_distances.txt", head=T, sep=",")
file$pdist = file$min_dist/file$exon_length

nrow(file[file$min_dist < 70 & file$min_dist > 3 & file$clinvar == 1,]) /nrow(file[file$clinvar == 1,])
nrow(file[file$min_dist < 70 & file$min_dist > 3 & file$clinvar == 0,]) /nrow(file[file$clinvar == 0,])


wilcox.test(file$min_dist[file$clinvar == 1], file$min_dist[file$X1000_genomes == 1])

clinvar <- file[file$clinvar == 1,]
kgenomes <- file[file$X1000_genomes == 1,]

# clinvar <- clinvar[clinvar$min_dist < quantile(clinvar$min_dist, 0.8), ]
# kgenomes <- kgenomes[kgenomes$min_dist < quantile(kgenomes$min_dist, 0.8), ]

clinvar = clinvar[clinvar$min_dist < 200,]
kgenomes = kgenomes[kgenomes$min_dist < 200,]

wilcox.test(clinvar$min_dist, kgenomes$min_dist)
t.test(clinvar$min_dist, kgenomes$min_dist)




compare_zs <- function(file1, file2, output_file) {
  file <- read.csv(file1, head=T)
  real = file[file$sim_id == "real",]
  sim = file[file$sim_id != "real",]
  
  call_3_z <- (real$all_0.3.bp - mean(sim$all_0.3.bp)) / sd(sim$all_0.3.bp)
  call_4_69_z <- (real$all_4.69.bp - mean(sim$all_4.69.bp)) / sd(sim$all_4.69.bp)
  call_70_z <- (real$all_70..bp - mean(sim$all_70..bp)) / sd(sim$all_70..bp)
  
  cexon_5_3_z <- (real$X5_0.3.bp - mean(sim$X5_0.3.bp)) / sd(sim$X5_0.3.bp)
  cexon_5_4_69_z <- (real$X5_4.69.bp - mean(sim$X5_4.69.bp)) / sd(sim$X5_4.69.bp)
  cexon_5_70_z <- (real$X5_70..bp - mean(sim$X5_70..bp)) / sd(sim$X5_70..bp)
  
  cexon_3_3_z <- (real$X3_0.3.bp - mean(sim$X3_0.3.bp)) / sd(sim$X3_0.3.bp)
  cexon_3_4_69_z <- (real$X3_4.69.bp - mean(sim$X3_4.69.bp)) / sd(sim$X3_4.69.bp)
  cexon_3_70_z <- (real$X3_70..bp - mean(sim$X3_70..bp)) / sd(sim$X3_70..bp)
  
  
  
  file <- read.csv(file2, head=T)
  real = file[file$sim_id == "real",]
  sim = file[file$sim_id != "real",]
  
  nall_3_z <- (real$all_0.3.bp - mean(sim$all_0.3.bp)) / sd(sim$all_0.3.bp)
  nall_4_69_z <- (real$all_4.69.bp - mean(sim$all_4.69.bp)) / sd(sim$all_4.69.bp)
  nall_70_z <- (real$all_70..bp - mean(sim$all_70..bp)) / sd(sim$all_70..bp)
  
  nexon_5_3_z <- (real$X5_0.3.bp - mean(sim$X5_0.3.bp)) / sd(sim$X5_0.3.bp)
  nexon_5_4_69_z <- (real$X5_4.69.bp - mean(sim$X5_4.69.bp)) / sd(sim$X5_4.69.bp)
  nexon_5_70_z <- (real$X5_70..bp - mean(sim$X5_70..bp)) / sd(sim$X5_70..bp)
  
  nexon_3_3_z <- (real$X3_0.3.bp - mean(sim$X3_0.3.bp)) / sd(sim$X3_0.3.bp)
  nexon_3_4_69_z <- (real$X3_4.69.bp - mean(sim$X3_4.69.bp)) / sd(sim$X3_4.69.bp)
  nexon_3_70_z <- (real$X3_70..bp - mean(sim$X3_70..bp)) / sd(sim$X3_70..bp)
  
  all_zs <- c(call_3_z, nall_3_z, call_4_69_z, nall_4_69_z, call_70_z, nall_70_z)
  start_zs <- c(cexon_5_3_z, nexon_5_3_z, cexon_5_4_69_z, nexon_5_4_69_z, cexon_5_70_z, nexon_5_70_z)
  end_zs <- c(cexon_3_3_z, nexon_3_3_z, cexon_3_4_69_z, nexon_3_4_69_z, cexon_3_70_z, nexon_3_70_z)
  
  disease_status <- c("ClinVar", "1000 Genomes", "ClinVar", "1000 Genomes", "ClinVar", "1000 Genomes")
  group <- c("0-3","0-3", "4-69", "4-69", "70+", "70+")
  
  library(ggplot2)
  library(reshape2)
  all <- data.frame(all_zs, disease_status, group)
  start <- data.frame(start_zs, disease_status, group)
  end <- data.frame(end_zs, disease_status, group)
  
  data_all <- melt(all, id=c("group", "disease_status"))
  data_start <- melt(start_zs, id=c("group", "disease_status"))
  data_end <- melt(end_zs, id=c("group", "disease_status"))
  
  
  upper = ceiling(max(3.5, max(call_3_z, call_4_69_z, call_70_z, cexon_3_3_z, cexon_3_4_69_z, cexon_3_70_z, cexon_5_3_z, cexon_5_4_69_z, cexon_5_70_z, nall_3_z, nall_4_69_z, nall_70_z, nexon_3_3_z, nexon_3_4_69_z, nexon_3_70_z, nexon_5_3_z, nexon_5_4_69_z, nexon_5_70_z)))
  lower = floor(min(-3.5, min(call_3_z, call_4_69_z, call_70_z, cexon_3_3_z, cexon_3_4_69_z, cexon_3_70_z, cexon_5_3_z, cexon_5_4_69_z, cexon_5_70_z, nall_3_z, nall_4_69_z, nall_70_z, nexon_3_3_z, nexon_3_4_69_z, nexon_3_70_z, nexon_5_3_z, nexon_5_4_69_z, nexon_5_70_z)))
  
  
  plot_zs <- function(data, upper, lower, main) {
    
    
    theme_Publication <- function(base_size=14) {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size)
        + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text.x = element_text(size=rel(1)),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                # legend.margin = unit(0, "cm"),
                legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    }
    
    plot <- ggplot(data, aes(x=group, fill=disease_status, y=value)) +   
      geom_bar(stat="identity", position=position_dodge()) + 
      labs(x="Region", y="Z", title=main) +
      scale_fill_manual(values = c("red", "blue")) + 
      scale_y_continuous(limits=c(lower, upper), breaks=seq(lower, upper, 1)) + 
      geom_hline(yintercept = 1.96, lty=2) +
      geom_hline(yintercept = -1.96, lty=2) +
      geom_hline(yintercept = 0, lty=1) +
      theme_Publication() + 
      theme(legend.title=element_blank())
    return(plot)
  }
  
  
  all <- plot_zs(data_all, upper, lower, "Full Exon")
  start <- plot_zs(data_start, upper, lower, "5'")
  end <- plot_zs(data_end, upper, lower, "3'")
  
  library(gridExtra)
  library(grid)
  plot <- grid.arrange(all, start, end)
  ggsave(output_file, plot, width=7, height=12)
}

file1 <- "results/clinvar/clinvar_simulations.csv"
file2 <- "results/clinvar/1000_genomes_simulations.csv"
output <- "results/clinvar/simulations.pdf"
compare_zs(file1, file2, output)

file1 <- "results/clinvar/clinvar_simulations_ese_overlaps.csv"
file2 <- "results/clinvar/1000_genomes_simulations_ese_overlaps.csv"
output <- "results/clinvar/simulations_ese_overlaps.jpg"
compare_zs(file1, file2, output)


# plot the z scores for the number of hits on an ese compared to randomly simulated
# positions in the 4-69 bp window at exon ends
clinvar_ese_hit <- read.csv("results/clinvar/clinvar_ese_hit_simulation_4_69.csv", head=T)
real <- clinvar_ese_hit[clinvar_ese_hit$simulation == "real",]
sims <- clinvar_ese_hit[clinvar_ese_hit$simulation != "real",]
cz5 <- (real$X5_ese_hit_count - mean(sims$X5_ese_hit_count)) / sd(sims$X5_ese_hit_count)
cz3 <- (real$X3_ese_hit_count - mean(sims$X3_ese_hit_count)) / sd(sims$X3_ese_hit_count)
clinvar_ese_hit <- read.csv("results/clinvar/1000_genomes_ese_hit_simulation_4_69.csv", head=T)
real <- clinvar_ese_hit[clinvar_ese_hit$simulation == "real",]
sims <- clinvar_ese_hit[clinvar_ese_hit$simulation != "real",]
gz5 <- (real$X5_ese_hit_count - mean(sims$X5_ese_hit_count)) / sd(sims$X5_ese_hit_count)
gz3 <- (real$X3_ese_hit_count - mean(sims$X3_ese_hit_count)) / sd(sims$X3_ese_hit_count)

disease_status <- c("ClinVar", "1000 Genomes", "ClinVar", "1000 Genomes")
group <- c("5'","5'", "3'", "3'")
zs <- c(cz5, gz5, cz3, gz3)
all <- data.frame(zs, disease_status, group)
pvals <- paste("p = ", signif(2*pnorm(-abs(zs)), digits=3), sep="")


library(reshape2)
data_all <- melt(all, id=c("group", "disease_status"))
data_all$group <- factor(data_all$group, levels = c("5'", "3'"))

lower = floor(min(zs))
upper = ceiling(max(zs))
ggplot(data_all, aes(x=group, fill=disease_status, y=value)) +   
  geom_bar(stat="identity", position=position_dodge()) + 
  labs(x="Exon end", y="Z", title="PTC ESE hits in 4-69 bp region") +
  scale_fill_manual(values = c("red", "blue")) + 
  scale_y_continuous(limits=c(lower, upper), breaks=seq(lower, upper, 1)) + 
  geom_hline(yintercept = 1.96, lty=2) +
  geom_hline(yintercept = 0, lty=1) +
  geom_text(aes(label=pvals), position=position_dodge(width=0.9), vjust=-0.45, cex=3) +
  theme_Publication() + 
  theme(legend.title=element_blank())

