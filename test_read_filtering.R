library(ggplot2)
library(reshape2)


file = "results/test_bam_filtering/final.csv"


data = read.csv(file, head = T)
head(data)

wilcox.test(data$full_prop, data$intersect_prop, paired = T)

data.subset = data.frame(data$sample, data$full_prop, data$intersect_prop)

melt.data = melt(data.subset)
head(melt.data)

geom_bar(data = melt.data, aes())

plot = ggplot(data = melt.data, aes(x = data.sample, y = value, fill = variable)) +
  geom_col(position = "dodge") +
  # geom_vline(xintercept = 1.5, lty = 2, col = "#aaaaaa") +
  # geom_vline(xintercept = 2.5, lty = 2, col = "#aaaaaa") +
  scale_fill_manual(values = c("RoyalBlue", "#e74b4f"), labels = c("Full", "Exon junction intersect")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Sample ID", y = "Proportion of reads retained") +
  theme_minimal() +
  theme(
    legend.title = element_blank()
  )

ggsave(plot, filename = "results/graphs/bam_filtering_comparison.pdf", width = 8, height = 5)
