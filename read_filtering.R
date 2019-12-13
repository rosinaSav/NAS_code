library(ggplot2)
library(reshape2)


file = "results/bam_filtering/final.csv"


data = read.csv(file, head = T)
head(data)

wilcox.test(data$full_prop, data$intersect_prop, paired = T)

data.subset = data.frame(data$sample, data$full_prop, data$intersect_prop)

melt.data = melt(data.subset)

plot = ggplot(data = melt.data, aes(x = data.sample, y = value, fill = variable)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#2c3e50", "#e67e22"), labels = c("Full", "Exon junction intersect")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Sample ID", y = "Proportion of reads retained\nafter filtering", fill = "Sample Type") +
  theme_minimal()

plot
ggsave(plot, filename = "results/graphs/bam_filtering_comparison.eps", width = 8, height = 5)
