install.packages("fitdistrplus")
library(fitdistrplus)
library(ggplot2)
library(MASS)
library(scales)


#### Analysis for random spacers ####

randomdata <- read.delim('../pipeline/spacersimulation_v2/random_set.txt', header = F, sep=' ')
randomdata <- t(randomdata)

comparedSpacers <- 430 #number of spacers simulated in each run of the simulation
simulationCycles <- nrow(randomdata) #number of simulations

randomdataFrame <- data.frame("clusters" = comparedSpacers - randomdata) #get number of clusters for each run
randomdataFrame$identity <- randomdataFrame$clusters/comparedSpacers #identity is calculated by dividing the number of clusters by the number of spacers compared
randomdataFrame$clusters <- randomdata

meanSim <- mean(randomdataFrame$identity)

##### Analysis for PAM spacers ####

PAMdata <- read.delim('../pipeline/spacersimulation_v2/PAM_set.txt', header = F, sep=' ')
PAMdata <- t(PAMdata)

PAMdataFrame <- data.frame("clusters" = comparedSpacers - PAMdata) #get number of clusters for each run
PAMdataFrame$identity <- PAMdataFrame$clusters/comparedSpacers #identity is calculated by dividing the number of clusters by the number of spacers compared
PAMdataFrame$clusters <- PAMdata

meanSimPAM <- mean(PAMdataFrame$identity)

# Fit data into normal dsitribution
fitRandom <- fitdistr(randomdataFrame$identity, densfun="normal")    #Fit the simulated data to a normal distribution
fitPAM <- fitdistr(PAMdataFrame$identity, densfun="normal")

#Plot 
ybreaks = seq(0,50,5) 

n = simulationCycles
meanRandom = fitRandom$estimate[1]
sdRandom = fitRandom$estimate[2]
meanPAM = fitPAM$estimate[1]
sdPAM = fitPAM$estimate[2]
observedPooled <- 0.58
bw=0.001
binwidth = bw # passed to geom_histogram and stat_function
bins = 100
randomdataTheoretical <- data.frame(x = rnorm(simulationCycles, fitRandom$estimate[1], fitRandom$estimate[2]))
randomdataPAMTheoretical <- data.frame(y = rnorm(simulationCycles, fitPAM$estimate[1], fitPAM$estimate[2]))
observedTypeVI <- 0.58

fontSize = 3
fontFamily = "Arial"
colorpalette <- c("#E69F00", "#999999", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(randomdataFrame, aes(x=identity)) + 
  geom_histogram(aes(y=..density..), bins=bins, fill=colorpalette[1], alpha = 0.8) +
  geom_histogram(data = PAMdataFrame, aes(y=..density..), bins=bins, fill=colorpalette[2], alpha=0.8) +
  scale_colour_manual(values=colorpalette) +
  xlim(0,1) +
  theme_minimal() +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0), breaks = pretty_breaks(n=10)) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_function(fun = function(x) dnorm(x, mean = fitRandom$estimate[1], sd = fitRandom$estimate[2]),
                color = colorpalette[1], size = 0.6, alpha = 0.6, n = 1000) +
  stat_function(fun = function(x) dnorm(x, mean = fitPAM$estimate[1], sd = fitPAM$estimate[2]),
                color = colorpalette[2], size = 0.6, alpha = 0.6, n = 1000) +
  annotate("text", label = "Type VI\nspacer pool", x = observedTypeVI, y=9, size = fontSize, colour = "black", family = fontFamily)  +
  # annotate("text", label = "Type VI\n(rep. b)", x = observedB+0.045, y=38, size = fontSize, colour = "black", family = fontFamily)  +
  #  annotate("text", label = "Type VI\n (rep. e)", x = observedE-0.045, y=38, size = fontSize, colour = "black", family = fontFamily)  +
  annotate("text", label = "Random\nspacers", x = fitRandom$estimate[1]+0.05, y=20, size = fontSize, colour = colorpalette[1], family = fontFamily)  +
  annotate("text", label = "Random\nPAM spacers", x = fitPAM$estimate[1], y=28, size = fontSize, colour = colorpalette[2], family = fontFamily)  +
  #geom_point(stat="identity") +
  #annotate("point", x=0.6, y=0, label="b", color="red")+
  annotate("pointrange", x=observedTypeVI, y=0, ymin=0, ymax=4, color="black")+
  #annotate("pointrange", x=observedB, y=0, ymin=0, ymax=50, color="black")+
  #annotate("pointrange", x=observedE, y=0, ymin=0, ymax=50, color="black")+
  #annotate(x=0.6, ymin=0, ymax=0, color="red") +
  #geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))
  xlab("Proportion of identical spacers with the type II-C spacer pool") +
  theme(text=element_text(family=fontFamily, size = fontSize),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(size = fontSize*3),
        axis.text.x = element_text(size = fontSize*3),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.1,1,0.1,1),"cm")) #top, right, bottom, left

#### Statistical testing ####

P_value_randomSpacers <- dnorm(observedPooled, fitRandom$estimate[1], fitRandom$estimate[2]) #T-test to compare our observed value to theoretical distribution modeled on the simulated data
P_value_randomSpacers
P_value_PAMSpacers <- dnorm(observedPooled, fitPAM$estimate[1], fitPAM$estimate[2]) #T-test to compare our observed value to theoretical distribution modeled on the simulated data
P_value_PAMSpacers
