file1 <- read.table("./results/disease_snps/disease__analysis_ptc_info.txt", head=T, sep=",")
hist(file1$min_dist, breaks=seq(min(file1$min_dist),max(file1$min_dist)+10, by=10), freq=F)

file2 <- read.table("./results/disease_snps/disease__analysis_other_info.txt", head=T, sep=",")
hist(file2$min_dist, breaks=seq(min(file2$min_dist),max(file2$min_dist)+10, by=10), freq = F)

nrow(file1[file1$min_dist<70,])/nrow(file1)
nrow(file2[file2$min_dist<70,])/nrow(file2)
