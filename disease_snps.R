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

read_in_simulations = function(simulant_prefix, simulant_number, ids, column_names) {
  sim_mat = matrix(NA, 1, 19)
  # colnames(sim_mat) = column_names
  out_list = rep(list(sim_mat), length(ids))
  for (sim in 1:simulant_number) {
    print(sim)
    current_data = read.table(paste(simulant_prefix, sim, ".bed", sep = ""), sep="\t"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    for (exon in 1:dim(current_data)[1]) {
      current_id = current_data[exon, "id"]
      if (current_id %in% ids) {
        out_list[[match(current_id, ids)]] = rbind(out_list[[match(current_id, ids)]], current_data[exon, ])
      }
    }
  }
  return(out_list)
}