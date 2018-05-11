file <- read.table("./results/disease_snps/disease__analysis_ptc_info.txt", head=T, sep=",")
hist(file$min_dist, breaks=seq(min(file$min_dist),max(file$min_dist)+10, by=10), freq=F)

file <- read.table("./results/disease_snps/disease__analysis_other_info.txt", head=T, sep=",")
hist(file$min_dist, breaks=seq(min(file$min_dist),max(file$min_dist)+10, by=10), freq = F)

nrow(file[file$min_dist<70,])/nrow(file)
nrow(file[file$min_dist<70,])
