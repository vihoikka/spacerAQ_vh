install.packages("ggExtra")
install.packages("Hmisc")
install.packages("tolerance-package")
install.packages("epitools")
install.packages("Rmisc")

isOnORF <<- function(spacerlocation,strandedness,genomeFunc) {
  strandedness = as.character(strandedness)
  genome <- subset(genomeFunc, dir != strandedness)
  genome <- subset(genome, min <= (spacerlocation))
  genome <- subset(genome, min >= (spacerlocation - max(genome[,5])))
  if (nrow(genome) == 0) {return (FALSE) }
  for (i in 1:nrow(genome)) {
    orfStart <- (genome[i,3]) #Only then check positions
    orfStop <- (genome[i,4])
    if (isTRUE(spacerlocation >= orfStart & isTRUE(spacerlocation < (orfStop-30)))) { #If the spacer is on an ORF, check if their direction is the same
      return (TRUE)
    }
  }
  return (FALSE)
}

library(ggplot2)
library(Hmisc)
library(ggExtra)
library(grid)
library(gridExtra)
library(ggrepel)
library(epitools)
library(Rmisc)

dataReplicates <- read.csv('../pipeline/ORFs/ORF_spacer_stats.csv', header = T, sep=',')

dataReplicates$totalExcludeIntergenics <- dataReplicates[,"total_hits"]-dataReplicates[,"intergenics"] #We only count intergenic spacers


dataReplicates$r_ratio <- dataReplicates[,"mRNA_hits"]/dataReplicates[,"totalExcludeIntergenics"]
dataReplicates$nonmRNA_hits <- dataReplicates[,"totalExcludeIntergenics"]-dataReplicates[,"mRNA_hits"]

dataReplicates$mRNAtoNonmRNAratio <- dataReplicates[,"mRNA_hits"]/dataReplicates[,"nonmRNA_hits"]

dataReplicatesPhage <- subset(dataReplicates, target=="fcl2")
dataReplicatesSelf <- subset(dataReplicates, target=="b185")

dataReplicatesC1 <- subset(dataReplicates, locus=="C1")
dataReplicatesC2 <- subset(dataReplicates, locus=="C2")

dataReplicatesList <- list(dataReplicatesPhage,dataReplicatesSelf)


# Fetch pooled data
dataPooled <- read.csv('../pipeline/Pooled_spacers/ORF_targets/ORF_spacer_stats_pooled.csv', header = T, sep=',')
dataPooled$totalExcludeIntergenics <- dataPooled[,"total_hits"]-dataPooled[,"intergenics"]
dataPooled$r_ratio <- dataPooled[,"mRNA_hits"]/dataPooled[,"totalExcludeIntergenics"]
dataPooled$nonmRNA_hits <- dataPooled[,"totalExcludeIntergenics"]-dataPooled[,"mRNA_hits"]
dataPooled$mRNAtoNonmRNAratio <- dataPooled[,"mRNA_hits"]/dataPooled[,"nonmRNA_hits"]

dataPooledC1 <- subset(dataPooled, locus=="C1")
dataPooledC2 <- subset(dataPooled, locus=="C2")

genomes <- c("../pipeline/FCL2_orfs.csv","../pipeline/b185_orfs.csv")
genomeLengths <- c(47142,3261404)
genomeCounter <- 1 #For determining what to plot
Z <- 1000 #number of runs

ORFStatsList<-list(fcl2=vector(mode="numeric", length=Z), b185=vector(mode="numeric", length=Z))

##### Stats

binomValues_C1_phage <- seq(0,dataPooledC1$totalExcludeIntergenics[1],by =1)
C1_phage_function <- dbinom(binomValues_C1_phage, dataPooledC1$totalExcludeIntergenics[1], 0.5)
C1_phage_probFrame <- data.frame(x=binomValues_C1_phage, y=C1_phage_function)
C1_phage_dens <- data.frame(y=C1_phage_probFrame$y/sum(C1_phage_probFrame$y)*100, x=seq(0,1, by=1/dataPooledC1$totalExcludeIntergenics[1]))

binomValues_C1_self <- seq(0,dataPooledC1$totalExcludeIntergenics[2],by =1)
C1_self_function <- dbinom(binomValues_C1_self, dataPooledC1$totalExcludeIntergenics[2], 0.5)
C1_self_probFrame <- data.frame(x=binomValues_C1_self, y=C1_self_function)
C1_self_dens <- data.frame(y=C1_self_probFrame$y/sum(C1_self_probFrame$y)*100, x=seq(0,1, by=1/dataPooledC1$totalExcludeIntergenics[2]))

binomValues_C2_phage <- seq(0,dataPooledC2$totalExcludeIntergenics[1],by =1)
C2_phage_function <- dbinom(binomValues_C2_phage, dataPooledC2$totalExcludeIntergenics[1], 0.5)
C2_phage_probFrame <- data.frame(x=binomValues_C2_phage, y=C2_phage_function)
C2_phage_dens <- data.frame(y=C2_phage_probFrame$y/sum(C2_phage_probFrame$y)*100, x=seq(0,1, by=1/dataPooledC2$totalExcludeIntergenics[1]))

binomValues_C2_self <- seq(0,dataPooledC2$totalExcludeIntergenics[2],by =1)
C2_self_function <- dbinom(binomValues_C2_self, dataPooledC2$totalExcludeIntergenics[2], 0.5)
C2_self_probFrame <- data.frame(x=binomValues_C2_self, y=C2_self_function)
C2_self_dens <- data.frame(y=C2_self_probFrame$y/sum(C2_self_probFrame$y)*100, x=seq(0,1, by=1/dataPooledC2$totalExcludeIntergenics[2]))

# Confidence intervals

C1_phage_95conf <- binconf(x=0.5*dataPooledC1$totalExcludeIntergenics[1], n=dataPooledC1$totalExcludeIntergenics[1], alpha = 0.025)
C1_self_95conf <- binconf(x=0.5*dataPooledC1$totalExcludeIntergenics[2], n=dataPooledC1$totalExcludeIntergenics[2])
C2_phage_95conf <- binconf(x=0.5*dataPooledC2$totalExcludeIntergenics[1], n=dataPooledC2$totalExcludeIntergenics[1])
C2_self_95conf <- binconf(x=0.5*dataPooledC2$totalExcludeIntergenics[2], n=dataPooledC2$totalExcludeIntergenics[2])

C1_phageProbablity <- pbinom(subset(dataPooledC1, target == "fcl2")$r_ratio*dataPooledC1$totalExcludeIntergenics[1], dataPooledC1$totalExcludeIntergenics[1],0.5)
C1_selfProbablity <- pbinom(subset(dataPooledC1, target == "b185")$r_ratio*dataPooledC1$totalExcludeIntergenics[2], dataPooledC1$totalExcludeIntergenics[2],0.5)
C2_phageProbablity <- pbinom(subset(dataPooledC2, target == "fcl2")$r_ratio*dataPooledC2$totalExcludeIntergenics[1], dataPooledC2$totalExcludeIntergenics[1],0.5)
C2_selfProbablity <- pbinom(subset(dataPooledC2, target == "b185")$r_ratio*dataPooledC2$totalExcludeIntergenics[2], dataPooledC2$totalExcludeIntergenics[2],0.5)

# P-values

PvalueTable <- array(c(C1_selfProbablity,C1_phageProbablity,C2_selfProbablity,C2_phageProbablity),dim = c(2,2),dimnames = list(c("Self","Phage"),c("II-C","VI-B")))
PvalueTable

##### Common plot

mean_hist_fcl2 <- mean(ORFStatsList[[1]]) #Get mean of simulated data
mean_hist_b185 <- mean(ORFStatsList[[2]]) #Get mean of simulated data

SD_fcl2 <- sd(ORFStatsList[[1]])
SD_b185 <- sd(ORFStatsList[[2]])

combined <- rbind(dataReplicatesList[[1]],dataReplicatesList[[2]])

colorpalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

b185_color = colorpalette[1]
fcl2_color = colorpalette[2]


levels(combined$target)[levels(combined$target)=="b185"] <- "Self"
levels(combined$target)[levels(combined$target)=="fcl2"] <- "Phage"

#Pooled values
C1_phage_prop <- subset(dataPooledC1, target == "fcl2")$r_ratio
C1_self_prop <- subset(dataPooledC1, target == "b185")$r_ratio
C1_props <- c(C1_phage_prop,C1_self_prop)
C2_phage_prop <- subset(dataPooledC2, target == "fcl2")$r_ratio
C2_self_prop <- subset(dataPooledC2, target == "b185")$r_ratio


#II-C plot
C1plot <- ggplot(dataReplicatesC1, aes(x=target, y=r_ratio)) +
  scale_color_manual(values = c("black","black")) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_fill_manual(values = c(b185_color, fcl2_color)) +
  geom_point(data = dataPooledC1, aes(y=r_ratio, x = target, fill = target), size = 3, stroke = 1, shape = 21, position = position_dodge(width = 0.3), alpha = 1) +
  geom_point(size = 1, colour="black", stroke = 1, shape = 20, alpha = 0.8) +
  annotate("rect", xmin = 0, xmax = 3, ymin = C1_phage_95conf[2], ymax = C1_phage_95conf[3],
           alpha = .2, fill = fcl2_color)  +
  annotate("rect", xmin = 0, xmax = 3, ymin = C1_self_95conf[2], ymax = C1_self_95conf[3],
           alpha = .2, fill = b185_color)  +
  ggtitle("Phage mRNA targeting ratio") + 
  geom_hline(yintercept = mean(C1_phage_dens$x), linetype="dotted", color = fcl2_color, size = 0.8, alpha = 0.6) +
 # geom_hline(yintercept = mean_hist_b185, linetype="dotted", color = b185_color, size = 1.1, alpha = 0.6) +
  labs(y = "Proportion of mRNA-targeting spacers", x = "II-C", size = 16) +
  #scale_x_continuous(expand = c(0, 0)) +
  labs(tag = "A") +
  scale_x_discrete(
    labels=c("Self", "Phage")) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
       # legend.position = c(0,1),
        legend.justification = c(0, 1),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
       # axis.title.x = element_text("II-C spacers"),
        #plot.title = element_text(size= 11, hjust = 0.5),
        plot.title = element_blank(),
        plot.margin=unit(c(0.5,-0.2,0.5,0.5),"cm")) #top, right, bottom, left))

#VI-B plot
C2plot <- ggplot(dataReplicatesC2, aes(x=target, y=r_ratio)) +
  scale_color_manual(values = c("black","black")) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_fill_manual(values = c(b185_color, fcl2_color)) +
  geom_point(data = dataPooledC2, aes(y=r_ratio, x = target, fill = target), size = 3, stroke = 1, shape = 21, position = position_dodge(width = 0.3), alpha = 1) +
  geom_point(size = 1, colour="black", stroke = 1, shape = 20, alpha = 0.8) +
  ggtitle("Phage mRNA targeting ratio") + 
  annotate("rect", xmin = 0, xmax = 3, ymin = C2_self_95conf[2], ymax = C2_self_95conf[3],
           alpha = .2, fill = b185_color)  +  
  annotate("rect", xmin = 0, xmax = 3, ymin = C2_phage_95conf[2], ymax = C2_phage_95conf[3],
           alpha = .2, fill = fcl2_color)  +
  geom_hline(yintercept = mean(C1_phage_dens$x), linetype="dotted", color = fcl2_color, size = 0.8, alpha = 0.6) +
#  geom_hline(yintercept = mean_hist_fcl2, linetype="dotted", color = fcl2_color, size = 1.1, alpha = 0.6) +
 # geom_hline(yintercept = mean_hist_b185, linetype="dotted", color = b185_color, size = 1.1, alpha = 0.6) +
  labs(y = "", x = "VI-B", size = 16) +
  labs(tag = "B") +
  #scale_x_continuous(expand = c(0, 0)) +
  scale_x_discrete(
    labels=c("Self", "Phage")) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
       #legend.position = c(0,1),
        legend.justification = c(0, 1),
        legend.text = element_text(size = 12),
        axis.title.y = element_blank(),
       axis.text.y = element_text(size = 12),
       # axis.ticks = element_blank(),
       #axis.title.x = element_text("VI-B spacers"),
        #plot.title = element_text(size= 11, hjust = 0.5),
        plot.title = element_blank(),
        plot.margin=unit(c(0.5,-0.2,0.5,0.5),"cm")) #top, right, bottom, left))


hist_C1 <- ggplot(data.frame(x=c(0, 2)), aes(x)) +
  scale_colour_manual(name="Colors",values=c(Self=b185_color, Phage=b185_color)) +
  #geom_histogram(aes(ORFStatsList[[1]], fill = "Phage"),binwidth = 0.01, alpha = 0.5, fill = fcl2_color)+
 # geom_histogram(aes(ORFStatsList[[2]], fill = "Self"),binwidth = 0.01, alpha = 0.5, fill = b185_color)+
 # stat_function(fun = C1_phage_function) +
  geom_area(data = C1_phage_probFrame, aes(x = x/dataPooledC1$totalExcludeIntergenics[1], y = y*dataPooledC1$totalExcludeIntergenics[1]), stat = "identity", alpha = 0.5, fill = fcl2_color, n=20000) +
  geom_area(data = C1_self_probFrame, aes(x = x/dataPooledC1$totalExcludeIntergenics[2], y = y*dataPooledC1$totalExcludeIntergenics[2]), stat = "identity", alpha = 0.5, fill = b185_color, n=20000) +
  coord_flip(xlim = c(0,1)) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme_minimal() +
  labs(y = " ", x = " ", size = 16, tag="") +
  #geom_vline(xintercept = mean_hist, linetype  = "dotted") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "#FFFFFF", size=12),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(colour = "#FFFFFF"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.margin=unit(c(0.5,0,0.5,0),"cm")) #top, right, bottom, left)

hist_C2 <- ggplot()+
  scale_colour_manual(name="Colors",values=c(Self=b185_color, Phage=b185_color)) +
 # geom_histogram(aes(ORFStatsList[[1]], fill = "Phage"),binwidth = 0.01, alpha = 0.5, fill = fcl2_color)+
 # geom_histogram(aes(ORFStatsList[[2]], fill = "Self"),binwidth = 0.01, alpha = 0.5, fill = b185_color)+
  geom_area(data = C2_phage_probFrame, aes(x = x/dataPooledC2$totalExcludeIntergenics[1], y = y*dataPooledC2$totalExcludeIntergenics[1]), stat = "identity", alpha = 0.5, fill = fcl2_color, n=20000) +
  geom_area(data = C2_self_probFrame, aes(x = x/dataPooledC2$totalExcludeIntergenics[2], y = y*dataPooledC2$totalExcludeIntergenics[2]), stat = "identity", alpha = 0.5, fill = b185_color, n=20000) +
  coord_flip(xlim = c(0,1)) + 
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme_minimal() +
  labs(y = "", x = "", size = 16, tag="") +
  #geom_vline(xintercept = mean_hist, linetype  = "dotted") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "#FFFFFF", size=12),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(colour = "#FFFFFF"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.margin=unit(c(0.5,0,0.5,0),"cm")) #top, right, bottom, left)

grid.arrange(C1plot, hist_C1, C2plot, hist_C2, nrow = 1, ncol = 4, widths=c(4,1,4,1))

