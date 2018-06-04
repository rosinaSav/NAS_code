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

