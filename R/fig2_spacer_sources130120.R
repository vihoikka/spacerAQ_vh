install.packages("ggrepel")
library(ggplot2)
library(gridExtra)
library(ggrepel)
install.packages("viridis")
library(viridis)
library(melt)
library(reshape2)

# This script takes the target file produced by the pipeline (targets.csv). After making it long format, it plots the number of spacers
#originating from self, phage and other.

#Read the targets file
dada <- read.csv('../pipeline/targets.csv', header = TRUE, sep = ",")
dada <- subset(dada, select = c(sample, locus, phage, self))

#Make it long
dada2long <- melt(dada)

#Some tuning of labels
C1 <- subset(dada, locus=="C1")
C2 <- subset(dada, locus=="C2")

C1 <- subset(dada2long, locus=="C1")
C2 <- subset(dada2long, locus=="C2")

s <- ggplot(C1, aes(x=target,y=amount))
p <- ggplot(C2, aes(x=target,y=amount))

C1$replicate <- c("b","e","b","e")
C2$replicate <- c("b","d","e","b","d","e")

#Set y limit for plot
ylimit <- 600

fontsize = 20
globalFont = "Arial"
colorpalette <- c("#E69F00", "#999999", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Create plots for both loci
C1plot <- ggplot(C1, aes(x=variable,y=value)) +
  geom_bar(data=C1, aes(x=variable,y=value, fill=variable), stat = "summary", fun.y = "mean", alpha = 1, inherit.aes = FALSE) +
  geom_dotplot(binaxis='y', binwidth = 4, dotsize = 4,labels=sample) + ylab("Unique spacers") + xlab("Type II-C") + 
  geom_label_repel(aes(label = replicate), size = fontsize/3, segment.alpha = 0.5, nudge_x = 0.3) +
  scale_fill_manual(values=colorpalette) +
  scale_y_continuous(expand = c(0, 0)) +
  #geom_line(aes(x=variable,y=value, group=sample), alpha = 0.15) +
  coord_cartesian(ylim = c(0,ylimit)) +
  theme(text=element_text(size=fontsize)) +
  #eliminates background, gridlines, and chart border
  theme_minimal() +
  theme(text=element_text(family=globalFont),
        axis.line = element_line(color = 'black'), 
        axis.text=element_text(size=fontsize, color = 'black'),
        axis.title = element_text(size = fontsize),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), #top, right, bottom, left
        panel.grid.minor.x = element_blank())

C2plot <- ggplot(C2, aes(x=variable,y=value)) +
  geom_bar(data=C2, aes(x=variable,y=value, fill=variable), stat = "summary", fun.y = "mean", alpha = 1, inherit.aes = FALSE) +
  geom_dotplot(binaxis='y', binwidth = 4, dotsize = 4,labels=sample) + ylab("Unique spacers") + xlab("Type VI-B") + 
  geom_label_repel(aes(label = replicate), size = fontsize/3, segment.alpha = 0.5, nudge_x = 0.3) +
  scale_fill_manual(values=colorpalette) +
  #geom_line(aes(x=variable,y=value, group=sample), alpha = 0.15) +
  coord_cartesian(ylim = c(0,ylimit)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text=element_text(size=fontsize)) +
  #eliminates background, gridlines, and chart border
  theme_minimal() +
  theme(text=element_text(family=globalFont),
        axis.line = element_line(color = 'black'), 
        axis.text=element_text(size=fontsize, color = 'black'),
        axis.title = element_text(size = fontsize),
        axis.title.y = element_text(size=1, color = 'white'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=1, color = 'white'),
        legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), #top, right, bottom, left
        panel.grid.minor.x = element_blank())

#Draw both plots
grid.arrange(C1plot, C2plot, nrow = 1)

