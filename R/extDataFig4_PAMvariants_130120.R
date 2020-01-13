library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)

PAMdada <- read.csv('../pipeline/PAMs/PAMs.csv', header = TRUE, sep = ",")

PAMdada$target <- as.character(PAMdada$target)
PAMdada$target[PAMdada$target == "fcl2"] <- "Phage"
PAMdada$target[PAMdada$target == "b185"] <- "Self"

PAMdada$other <- PAMdada[,8]-(rowSums(PAMdada[,4:7]))

for(col in names(PAMdada[,c(3:7,9)])[-1]) {
  PAMdada[paste0(col, "_pct")] = PAMdada[col] / PAMdada[,8]
}

PAMdadalong <- melt(PAMdada[,c(-4,-5,-6,-7,-8,-9)],id.vars=c("sample","locus","target"))

PAMdadalongrev <- PAMdadalong[order(PAMdadalong$variable,decreasing=T),]


grid.newpage()

plot <- ggplot(PAMdadalong, aes(x=sample, y=value, fill=variable)) +
  scale_x_discrete(limits=c("B+P_b_+1_C1",
                            "B+P_e_+1_C1",
                            "B+P_b_+1_C2",
                            "B+P_d_+1_C2",
                            "B+P_e_+1_C2"),
                   labels=c("II-C b", 
                            "II-C e",
                            "VI-B b",
                            "VI-B d",
                            "VI-B e")) +
  scale_fill_brewer(palette="Accent",
                    labels=c("NNNNNTAAA (5bp gap)","NNNNTAAA (4bp gap)","NNNNNNTAAA (6bp gap)","NNNNNNNTAAA (7bp gap)","Other"),
                    direction = -1, #reverse order of colors
                    guide = guide_legend(reverse = TRUE)) +
  geom_bar(stat="identity", position = position_fill(reverse = TRUE)) +
  facet_grid(~target) +
  ggtitle("PAM variants") +
  ylab("Proportion") +
  guides(colour = guide_legend(reverse=T)) + #reverse order of stacking
  theme(axis.text.x = element_text(size=11),
        legend.position = "right",
        plot.title = element_text(hjust=0.5),
        legend.title = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14))
plot