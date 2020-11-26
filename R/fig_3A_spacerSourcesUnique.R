install.packages("ggrepel")
library(ggplot2)
library(gridExtra)
library(ggrepel)
install.packages("viridis")
library(viridis)
library(reshape2)

# This script takes the target file produced by the pipeline (targets.csv). After making it long format, it plots the number of spacers
#originating from self, phage and other.

#Read the targets file
dada <- read.csv('../pipeline/targets.csv', header = TRUE, sep = ",")
dada <- subset(dada, select = c(sample, locus, phage, self))

#Make it long
dada2long <- melt(dada)
aggr <- aggregate(value ~ sample, data = dada2long, FUN = sum)
dada2long <- merge(dada2long, aggr, by = "sample", all.x=TRUE, all.y=TRUE)
dada2long$proportion <- dada2long$value.x/dada2long$value.y
dada2long$replicate <- substring(dada2long$sample, 5, 5)

#fontsize = 20
colorpalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


samplesNames <- c(
  `C1` = "II-C",
  `C2` = "VI-B"
)

uniquespacers <- ggplot(data = dada2long, aes(x = variable, y = proportion)) +
  geom_bar(stat = "summary", fun = "mean") +
  geom_point() +
  geom_line((aes(group = replicate, color = replicate)), alpha = 0.4) +
  geom_label_repel(aes(label = replicate, fill = replicate), segment.alpha = 0.5, nudge_x = 0.3) +
  scale_fill_manual(values=colorpalette) +
  scale_color_manual(values=colorpalette) +
  facet_grid(~locus,labeller = as_labeller(samplesNames)) +
  xlab("Target genome") +
  ylab("Unique spacer count (proportion)") +
  theme_bw() +
  theme(
    legend.position = "none"
  )

ggsave(uniquespacers, filename = "fig_3A.png", height = 3, width = 3)

