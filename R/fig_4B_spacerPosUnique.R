#For plotting spacer position data. Takes in csv files, in which each row represents a position (usually a binning window ID)
#Used for Figure 3B

install.packages("gridExtra")
install.packages("ggplot2")

library(plyr) #For renaming a column
library(gridExtra)
library(grid)
library(reshape2) #For wide -> long
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)


#Fetch data for self hits in both loci and all replicates
b_C1_self_f <- read.csv('../pipeline/Binned/B+P_b_+1_C1_map_b185_binned_C1_mappedF.csv', header = TRUE, sep = ",")
e_C1_self_f <- read.csv('../pipeline/Binned/B+P_e_+1_C1_map_b185_binned_C1_mappedF.csv', header = TRUE, sep = ",")

b_C1_self_r <- read.csv('../pipeline/Binned/B+P_b_+1_C1_map_b185_binned_C1_mappedR.csv', header = TRUE, sep = ",")
e_C1_self_r <- read.csv('../pipeline/Binned/B+P_e_+1_C1_map_b185_binned_C1_mappedR.csv', header = TRUE, sep = ",")

b_C2_self_f <- read.csv('../pipeline/Binned/B+P_b_+1_C2_map_b185_binned_C2_mappedF.csv', header = TRUE, sep = ",")
d_C2_self_f <- read.csv('../pipeline/Binned/B+P_d_+1_C2_map_b185_binned_C2_mappedF.csv', header = TRUE, sep = ",")
e_C2_self_f <- read.csv('../pipeline/Binned/B+P_e_+1_C2_map_b185_binned_C2_mappedF.csv', header = TRUE, sep = ",")

b_C2_self_r <- read.csv('../pipeline/Binned/B+P_b_+1_C2_map_b185_binned_C2_mappedR.csv', header = TRUE, sep = ",")
d_C2_self_r <- read.csv('../pipeline/Binned/B+P_d_+1_C2_map_b185_binned_C2_mappedR.csv', header = TRUE, sep = ",")
e_C2_self_r <- read.csv('../pipeline/Binned/B+P_e_+1_C2_map_b185_binned_C2_mappedR.csv', header = TRUE, sep = ",")

#Fetch data for phage hits in both loci and all replicates
b_C1_phage_f <- read.csv('../pipeline/Binned/B+P_b_+1_C1_map_fcl2_binned_C1_mappedF.csv', header = TRUE, sep = ",")
e_C1_phage_f <- read.csv('../pipeline/Binned/B+P_e_+1_C1_map_fcl2_binned_C1_mappedF.csv', header = TRUE, sep = ",")

b_C1_phage_r <- read.csv('../pipeline/Binned/B+P_b_+1_C1_map_fcl2_binned_C1_mappedR.csv', header = TRUE, sep = ",")
e_C1_phage_r <- read.csv('../pipeline/Binned/B+P_e_+1_C1_map_fcl2_binned_C1_mappedR.csv', header = TRUE, sep = ",")

b_C2_phage_f <- read.csv('../pipeline/Binned/B+P_b_+1_C2_map_fcl2_binned_C2_mappedF.csv', header = TRUE, sep = ",")
d_C2_phage_f <- read.csv('../pipeline/Binned/B+P_d_+1_C2_map_fcl2_binned_C2_mappedF.csv', header = TRUE, sep = ",")
e_C2_phage_f <- read.csv('../pipeline/Binned/B+P_e_+1_C2_map_fcl2_binned_C2_mappedF.csv', header = TRUE, sep = ",")

b_C2_phage_r <- read.csv('../pipeline/Binned/B+P_b_+1_C2_map_fcl2_binned_C2_mappedR.csv', header = TRUE, sep = ",")
d_C2_phage_r <- read.csv('../pipeline/Binned/B+P_d_+1_C2_map_fcl2_binned_C2_mappedR.csv', header = TRUE, sep = ",")
e_C2_phage_r <- read.csv('../pipeline/Binned/B+P_e_+1_C2_map_fcl2_binned_C2_mappedR.csv', header = TRUE, sep = ",")

#Bacterial genome POI locations
C1_repeatSpacerArea_range <- c(755586,757197) #C1 repeat spacer array start & stop coordinates on B185
C1_repeatSpacerArea_middle <- sum(755586,757197)/2 #C1 repeat spacer array start & stop coordinates on B185

C2_repeatSpacerArea_range <- c(1284454,1285538) #C2 repeat spacer array start & stop coordinates on B185
C2_repeatSpacerArea_middle <- sum(1284454,1285538)/2 #C2 repeat spacer array start & stop coordinates on B185

C3_cas9 <- c(1945415,1949794)

#C1 existing phage-targeting
C1_existing_spacers_phage_F <- c(18513,  #c1s5
                                 23720, #c1s13
                                 12966 #c1s24
)

C1_existing_spacers_phage_F_mut <- c(42621) #c1s20

C1_existing_spacers_phage_R <- c(40160 #c1s23
) #c2s5
C1_existing_spacers_phage_F <- melt(C1_existing_spacers_phage_F)
C1_existing_spacers_phage_F$y = 0
C1_existing_spacers_phage_R <- melt(C1_existing_spacers_phage_R)
C1_existing_spacers_phage_R$y = 0

C1_existing_spacers_phage_F_mut <- melt(C1_existing_spacers_phage_F_mut)
C1_existing_spacers_phage_F_mut$y = 0

#C2 existing self spacers

C2_existing_spacers_self_F <- c(751992,  #c2s1
                                990469,  #c2s12
                                2006885, #c2s4
                                2730639, #c2s7
                                2626284) #c2s9
C2_existing_spacers_self_R <- c(1885737, #c2s2
                                1936610) #c2s5
C2_existing_spacers_self_F <- melt(C2_existing_spacers_self_F)
C2_existing_spacers_self_F$y = 0
C2_existing_spacers_self_R <- melt(C2_existing_spacers_self_R)
C2_existing_spacers_self_R$y = 0

#C2 existing phage spacers
C2_existing_spacers_phage_R <- c(8756, #c2s13
                                 44382, #c2s14
                                 40092, #c2s15
                                 42804) #c2s16
C2_existing_spacers_phage_R <- melt(C2_existing_spacers_phage_R)
C2_existing_spacers_phage_R$y = 0


#Fetch PAM data for phage and self
PAMs_FCL2 <- read.csv('../pipeline/PAMs_FCL2.csv', header = TRUE, sep = ";")
PAMs_FCL2_F <- subset(PAMs_FCL2,Direction=="forward")
PAMs_FCL2_R <- subset(PAMs_FCL2,Direction=="reverse")

PAMs_b185 <- read.csv('../pipeline/PAMs_b185.csv', header = TRUE, sep = ";")
PAMs_b185_F <- subset(PAMs_b185,Direction=="forward")
PAMs_b185_R <- subset(PAMs_b185,Direction=="reverse")

#Get average counts for bins for smooth graph
#Self
b_C1_f_self_avg <- b_C1_self_f$target
b_C1_r_self_avg <- b_C1_self_r$target
b_C2_f_self_avg <- b_C2_self_f$target
b_C2_r_self_avg <- b_C2_self_r$target

d_C2_f_self_avg <- d_C2_self_f$target
d_C2_f_self_avg <- d_C2_self_f$target

e_C1_f_self_avg <- e_C1_self_f$target
e_C1_r_self_avg <- e_C1_self_r$target
e_C2_f_self_avg <- e_C2_self_f$target
e_C2_r_self_avg <- e_C2_self_r$target

#Phage
b_C1_f_phage_avg <- b_C1_phage_f$target
b_C1_r_phage_avg <- b_C1_phage_r$target
b_C2_f_phage_avg <- b_C2_phage_f$target
b_C2_r_phage_avg <- b_C2_phage_r$target

d_C2_f_phage_avg <- d_C2_phage_f$target
d_C2_f_phage_avg <- d_C2_phage_f$target

e_C1_f_phage_avg <- e_C1_phage_f$target
e_C1_r_phage_avg <- e_C1_phage_r$target
e_C2_f_phage_avg <- e_C2_phage_f$target
e_C2_r_phage_avg <- e_C2_phage_r$target

# Calculate average bin hits over replicates on both strands
averageFs_self_C1 <- data.frame(b_C1_self_f$target,e_C1_self_f$target) #Combine data from replicates
averageRs_self_C1 <- data.frame(b_C1_self_r$target,e_C1_self_r$target)

averageFs_phage_C1 <- data.frame(b_C1_phage_f$target,e_C1_phage_f$target) #Combine data from replicates
averageRs_phage_C1 <- data.frame(b_C1_phage_r$target,e_C1_phage_r$target)

averageFs_self_C1$avg <- rowMeans(averageFs_self_C1) #Average the replicates into a new column
averageFs_self_C1$Pos <- b_C1_self_f$location #Add location...
averageFs_self_C1$binSize <- b_C1_self_f$binSize #And binsize information from one of the replicates

averageRs_self_C1$avg <- rowMeans(averageRs_self_C1) #Repeat for all replicates
averageRs_self_C1$Pos <- b_C1_self_r$location
averageRs_self_C1$binSize <- b_C1_self_r$binSize

averageFs_phage_C1$avg <- rowMeans(averageFs_phage_C1)
averageFs_phage_C1$Pos <- b_C1_phage_f$location
averageFs_phage_C1$binSize <- b_C1_phage_f$binSize

averageRs_phage_C1$avg <- rowMeans(averageRs_phage_C1)
averageRs_phage_C1$Pos <- b_C1_phage_r$location
averageRs_phage_C1$binSize <- b_C1_phage_r$binSize

#Repeat for C2
averageFs_self_C2 <- data.frame(b_C2_self_f$target,d_C2_self_f$target,e_C2_self_f$target)
averageRs_self_C2 <- data.frame(b_C2_self_r$target,d_C2_self_r$target,e_C2_self_r$target)

averageFs_phage_C2 <- data.frame(b_C2_phage_f$target,e_C2_phage_f$target) #Combine data from replicates
averageRs_phage_C2 <- data.frame(b_C2_phage_r$target,e_C2_phage_r$target)

averageFs_self_C2$avg <- rowMeans(averageFs_self_C2) #For explanation, see C1 above
averageFs_self_C2$Pos <- b_C2_self_f$location
averageFs_self_C2$binSize <- b_C2_self_f$binSize

averageRs_self_C2$avg <- rowMeans(averageRs_self_C2)
averageRs_self_C2$Pos <- b_C2_self_r$location
averageRs_self_C2$binSize <- b_C2_self_r$binSize

averageFs_phage_C2$avg <- rowMeans(averageFs_phage_C2)
averageFs_phage_C2$Pos <- b_C2_phage_f$location
averageFs_phage_C2$binSize <- b_C2_phage_f$binSize

averageRs_phage_C2$avg <- rowMeans(averageRs_phage_C2)
averageRs_phage_C2$Pos <- b_C2_phage_r$location
averageRs_phage_C2$binSize <- b_C2_phage_r$binSize

#Transform to long
averageFs_self_C1_long <- melt(averageFs_self_C1, id.vars=c("Pos","binSize"))
averageRs_self_C1_long <- melt(averageRs_self_C1, id.vars=c("Pos","binSize"))
averageFs_self_C2_long <- melt(averageFs_self_C2, id.vars=c("Pos","binSize"))
averageRs_self_C2_long <- melt(averageRs_self_C2, id.vars=c("Pos","binSize"))

averageFs_phage_C1_long <- melt(averageFs_phage_C1, id.vars=c("Pos","binSize"))
averageRs_phage_C1_long <- melt(averageRs_phage_C1, id.vars=c("Pos","binSize"))
averageFs_phage_C2_long <- melt(averageFs_phage_C2, id.vars=c("Pos","binSize"))
averageRs_phage_C2_long <- melt(averageRs_phage_C2, id.vars=c("Pos","binSize"))

#### Parameters ####

transparency <- 1/3
transparency2 <- 1/8
PAMtransparency <- 1/10
bordercolor <- "grey"
col1 <- colorpalette[1]
col2 <- colorpalette[2]
col3 <- colorpalette[3]

colorpalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

smoother1 = 1/3
smoother2 = 1/5

#ylimit = range(dadaR1$target)
ylimit1 = range(0,6.5)
ylimit2 = range(0,13)


ylimitORF = range(0,0.016)
ylimitORFself = range(0,0.006)

topLabelX <- 7*10^3
#topLabelXPhage <- 4*10^5
topLabelXPhage <- 47142/2
topLabelXSelf <- 3261404/2
topLabelYSelf <- ylimit2[2]-1.5
topLabelYPhage <- ylimit1[2]-1.5

locusLabelSeparator <- 1.8e5

globalSizeFactor = 1
lineThickness = 0.7 * globalSizeFactor
fontSizeTitle = 6 * globalSizeFactor
annotationSize = 5 * globalSizeFactor
existingSpacerSize = 1.5 * globalSizeFactor
spacerStroke = 0.6
tickSize = 15 * globalSizeFactor
gapBetweenPhageSelf = 0.3

transparency <- 1/3
transparency2 <- 1/4
PAMtransparency <- 1/10
bordercolor <- "grey"
col1 <- colorpalette[1]
col2 <- colorpalette[2]
col3 <- colorpalette[3]
smoothTrans <- 0.5
annoTrans <- 0.5
spacerColor <- "#D16103"
globalFont = "Arial"
rightSide = 0.2
leftSide = 0.2
topPhage = 0
topSelf = 0
bottomPhage = 0.2
bottomSelf = topPhage

human_numbers <- function(x = NULL, smbl ="", signif = 1){
  humanity <- function(y){
    
    if (!is.na(y)){
      tn <- round(abs(y) / 1e12, signif)
      b <- round(abs(y) / 1e9, signif)
      m <- round(abs(y) / 1e6, signif)
      k <- round(abs(y) / 1e3, signif)
      
      if ( y >= 0 ){
        y_is_positive <- ""
      } else {
        y_is_positive <- "-"
      }
      
      if ( k < 1 ) {
        paste0( y_is_positive, smbl, round(abs(y), signif ))
      } else if ( m < 1){
        paste0 (y_is_positive, smbl,  k , "kb")
      } else if (b < 1){
        paste0 (y_is_positive, smbl, m ,"mb")
      }else if(tn < 1){
        paste0 (y_is_positive, smbl, b ,"bn")
      } else {
        paste0 (y_is_positive, smbl,  comma(tn), "tn")
      }
    } else if (is.na(y) | is.null(y)){
      "-"
    }
  }
  
  sapply(x,humanity)
}
human_num   <- function(x){human_numbers(x, smbl = "")} 


C2_phage_plot_f <- ggplot(data = averageRs_phage_C2_long, aes(x = Pos+binSize/2, y = value),stat="identity") + #Average data as the base data
  geom_bar(data = b_C2_phage_r, stat = "identity", #Add geom_bars for each replicate to represent their bins
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col1) +
  scale_fill_brewer() +
  geom_bar(data = d_C2_phage_r, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col2) +
  geom_bar(data = e_C2_phage_r, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col3) +
  geom_line(data = PAMs_FCL2_F, stat = "density", aes(x=Minimum, y=..density..*10^5), adjust = 1/8, alpha = PAMtransparency*5, colour = "red", inherit.aes = F) + #Draw PAM frequency as red line
  geom_point(data = C2_existing_spacers_phage_R, aes(x=value, y=y), color = spacerColor, shape = 5, size = existingSpacerSize, stroke = spacerStroke) +
  geom_line(stat = 'smooth', method = 'loess', span = smoother1, size = lineThickness, alpha = 0.5, color = "black",formula = y ~ x,level = 0.999, se = FALSE) + #Draw the average as a smooth line
  annotate("text", label = "VI-B, phage", x = 0, y = topLabelYPhage, size = fontSizeTitle, colour = "black", hjust = 0)  +
  #geom_rug(data = PAMs_FCL2_F, aes(x=Minimum, y=0), position = position_jitter(height = 0), inherit.aes = F, alpha = transparency) + #This could be a spacer rug
  ylab("F strand (%)") +
  theme_minimal() +
  scale_y_continuous(position = "right", limits=c(0,ylimit1[2])) + #Flips y-values to right side
  scale_x_continuous(labels = human_num) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

C2_phage_plot_r <- ggplot(data = averageFs_phage_C2_long, aes(x = Pos+binSize/2, y = value),stat="identity") + #Average data as the base data
  scale_fill_brewer(palette="Accent") +
  geom_bar(data = b_C2_phage_f, stat = "identity", #Add geom_bars Ror each replicate to represent their bins
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col1) +
  geom_bar(data = d_C2_phage_f, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col2) +
  geom_bar(data = e_C2_phage_f, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col3) +
  geom_line(data = PAMs_FCL2_R, stat = "density", aes(x=Minimum, y=..density..*10^5), adjust = 1/8, alpha = PAMtransparency*5, colour = "red", inherit.aes = F) + #Draw PAM Rrequency as red line
  geom_line(stat = 'smooth', method = 'loess', span = smoother1, size = lineThickness, alpha = 0.5, color = "black",formula = y ~ x,level = 0.999, se = FALSE) + #Draw the average as a smooth line
  #geom_rug(data = PAMs_FCL2_R, aes(x=Minimum, y=0), binwidth = 10, position = position_jitter(height = 0), sides="top", inherit.aes = F, alpha = transparency) + #This could be a spacer rug
  ylab("R strand (%)") +
  theme_minimal() +
  #scale_y_reverse(position = "right") + #Turns the reverse plot upside down
  scale_y_reverse(limits=c(ylimit1[2],0), position = "right") + #Turns the reverse plot upside down
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


#### C2 SELF PLOTS ####
C2_self_plot_f <- ggplot(data = averageRs_self_C2_long, aes(x = Pos+binSize/2, y = value),stat="identity") + #Average data as the base data
  geom_bar(data = b_C2_self_r, stat = "identity", #Add geom_bars for each replicate to represent their bins
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col1, orientation = "x") +
  scale_fill_brewer() +
  geom_bar(data = d_C2_self_r, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col2, orientation = "x") +
  geom_bar(data = e_C2_self_r, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col3, orientation = "x") +
  geom_line(data = PAMs_b185_F, stat = "density", aes(x=Minimum, y=..density..*10^7), adjust = 1/8, alpha = PAMtransparency*5, colour = "red", inherit.aes = F) + #Draw PAM frequency as red line
  #geom_rug(data = PAMs_b185_F, aes(x=Minimum, y=0), position = position_jitter(height = 0), inherit.aes = F, alpha = transparency) + #This could be a spacer rug
  geom_line(stat = 'smooth', method = 'loess', span = smoother1, size = lineThickness, alpha = 0.5, color = "black",formula = y ~ x,level = 0.999, se = FALSE) + #Draw the average as a smooth line
  #geom_vline(xintercept = C1_repeatSpacerArea_middle, linetype="dotted") +
  # annotate("text", label = "VI-B", x = C2_repeatSpacerArea_middle+locusLabelSeparator, y = ylimit2[2], size = annotationSize, colour = "black", alpha = annoTrans)  +
  # geom_vline(xintercept = C2_repeatSpacerArea_middle, linetype="dotted") +
  # annotate("text", label = "II-C", x = C1_repeatSpacerArea_middle+locusLabelSeparator, y = ylimit2[2], size = annotationSize, colour = "black", alpha = annoTrans)  +
  #geom_vline(xintercept = 1997400, linetype="dotted", colour = "blue") +
  # annotate("text", label = "oriC", x = 1998700+1.5e5, y = ylimit2[2], size = annotationSize, colour = "black", alpha = annoTrans)  +
  #geom_vline(xintercept = mean(C3_cas9), linetype="dotted") +
  #geom_vline(xintercept = 1982000, linetype="dotted") +
  annotate("text", label = " VI-B, self", x = 0, y = topLabelYSelf, size = fontSizeTitle, colour = "black", hjust = 0, font = globalFont)  +
  geom_point(data = C2_existing_spacers_self_R, aes(x=value, y=y), color = spacerColor, size = existingSpacerSize, stroke = spacerStroke) +
  #ylab("F hits (%)") +
  scale_y_continuous(position = "right") + #Turns the reverse plot upside down
  scale_x_continuous(labels = human_num) +
  theme_minimal() +
  coord_cartesian(ylim = ylimit2) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

C2_self_plot_r <- ggplot(data = averageFs_self_C2_long, aes(x = Pos+binSize/2, y = value),stat="identity") + #Average data as the base data
  scale_fill_brewer(palette="Accent") +
  geom_bar(data = b_C2_self_f, stat = "identity", #Add geom_bars Ror each replicate to represent their bins
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col1) +
  geom_bar(data = d_C2_self_f, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col2) +
  geom_bar(data = e_C2_self_f, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col3) +
  geom_line(data = PAMs_b185_R, stat = "density", aes(x=Minimum, y=..density..*10^7), adjust = 1/8, alpha = PAMtransparency*5, colour = "red", inherit.aes = F) + #Draw PAM Rrequency as red line
  #geom_rug(data = PAMs_b185_R, aes(x=Minimum, y=0), binwidth = 10, position = position_jitter(height = 0), sides="top", inherit.aes = F, alpha = transparency) + #This could be a spacer rug
  geom_point(data = C2_existing_spacers_self_F, aes(x=value, y=y), color = spacerColor, size = existingSpacerSize, stroke = spacerStroke) +
  geom_line(stat = 'smooth', method = 'loess', span = smoother1, size = lineThickness, alpha = 0.5, color = "black",formula = y ~ x,level = 0.999, se = FALSE) + #Draw the average as a smooth line
  ylab("R targets (%)") +
  theme_minimal() +
  scale_y_reverse(limits=c(ylimit2[2],0), position = "right") + #Turns the reverse plot upside down
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


#### C1 PHAGE PLOTS ####
C1_phage_plot_f <- ggplot(data = averageRs_phage_C1_long, aes(x = Pos+binSize/2, y = value),stat="identity") + #Average data as the base data
  #ggtitle(" Type II-C, phage") +
  geom_bar(data = b_C1_phage_r, stat = "identity", #Add geom_bars for each replicate to represent their bins
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col1) +
  scale_fill_brewer() +
  geom_bar(data = e_C1_phage_r, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col3) +
  geom_line(data = PAMs_FCL2_F, stat = "density", aes(x=Minimum, y=..density..*10^5), adjust = 1/8, alpha = PAMtransparency*5, colour = "red", inherit.aes = F) + #Draw PAM frequency as red line
  geom_point(data = C1_existing_spacers_phage_R, aes(x=value, y=y), color = spacerColor, fill = spacerColor, shape = 23, size = existingSpacerSize, stroke = spacerStroke) +
  #geom_rug(data = PAMs_FCL2_F, aes(x=Minimum, y=0), position = position_jitter(height = 0), inherit.aes = F, alpha = transparency) + #This could be a spacer rug
  geom_line(stat = 'smooth', method = 'loess', span = smoother1, size = lineThickness, alpha = 0.5, color = "black",formula = y ~ x,level = 0.999, se = FALSE) + #Draw the average as a smooth line
  annotate("text", label = " II-C, phage", x = 0, y = topLabelYPhage, size = fontSizeTitle, colour = "black", hjust = 0)  +
  ylab("F targets (%)") +
  theme_minimal() +
  coord_cartesian(ylim = ylimit1) +
  scale_x_continuous(labels = human_num) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

C1_phage_plot_r <- ggplot(data = averageFs_phage_C1_long, aes(x = Pos+binSize/2, y = value),stat="identity") + #Average data as the base data
  scale_fill_brewer(palette="Accent") +
  geom_bar(data = b_C1_phage_f, stat = "identity", #Add geom_bars Ror each replicate to represent their bins
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col1, orientation = "x") +
  geom_bar(data = e_C1_phage_f, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col3, orientation = "x") +
  geom_line(data = PAMs_FCL2_R, stat = "density", aes(x=Minimum, y=..density..*10^5), adjust = 1/8, alpha = PAMtransparency*5, colour = "red", inherit.aes = F) + #Draw PAM Rrequency as red line
  geom_point(data = C1_existing_spacers_phage_F, aes(x=value, y=y), color = spacerColor, fill = "red", size = existingSpacerSize) +
  geom_point(data = C1_existing_spacers_phage_F_mut, aes(x=value, y=y), color = spacerColor, shape = 1, size = existingSpacerSize, stroke = spacerStroke) +
  #geom_rug(data = PAMs_FCL2_R, aes(x=Minimum, y=0), binwidth = 10, position = position_jitter(height = 0), sides="top", inherit.aes = F, alpha = transparency) + #This could be a spacer rug
  geom_line(stat = 'smooth', method = 'loess', span = smoother1, size = lineThickness, alpha = 0.5, color = "black",formula = y ~ x,level = 0.999, se = FALSE) + #Draw the average as a smooth line
  ylab("R targets (%)") +
  theme_minimal() +
  scale_y_reverse(limits=c(ylimit1[2],0)) + #Turns the reverse plot upside down
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

####C1 SELF PLOTS ####
C1_self_plot_f <- ggplot(data = averageRs_self_C1_long, aes(x = Pos+binSize/2, y = value),stat="identity") + #Average data as the base data
  #ggtitle(" Type II-C, self") +
  geom_bar(data = b_C1_self_r, stat = "identity", #Add geom_bars for each replicate to represent their bins
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col1, orientation = "x") +
  geom_bar(data = e_C1_self_r, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,colour=bordercolor, fill=col3, orientation = "x") +
  geom_line(data = PAMs_b185_F, stat = "density", aes(x=Minimum, y=..density..*10^7), adjust = 1/8, alpha = PAMtransparency*5, colour = "red", inherit.aes = F) + #Draw PAM frequency as red line
  #geom_rug(data = PAMs_b185_F, aes(x=Minimum, y=0), position = position_jitter(height = 0), inherit.aes = F, alpha = transparency) + #This could be a spacer rug
  #geom_vline(xintercept = C1_repeatSpacerArea_middle, linetype="dotted") +
  # annotate("text", label = "II-C", x = C1_repeatSpacerArea_middle+locusLabelSeparator, y = ylimit2[2], size = annotationSize, colour = "black", alpha = annoTrans)  +
  #geom_vline(xintercept = C2_repeatSpacerArea_middle, linetype="dotted") +
  # annotate("text", label = "VI-B", x = C2_repeatSpacerArea_middle+locusLabelSeparator, y = ylimit2[2], size = annotationSize, colour = "black", alpha = annoTrans)  +
  # geom_vline(xintercept = 1997400, linetype="dotted", colour = "blue") +
  geom_line(stat = 'smooth', method = 'loess', span = smoother1, size = lineThickness, alpha = 0.5, color = "black",formula = y ~ x,level = 0.999, se = FALSE) + #Draw the average as a smooth line
  # annotate("text", label = "oriC", x = 1998700+1.5e5, y = ylimit2[2], size = annotationSize, colour = "black", alpha = annoTrans)  +
  #geom_point(data = C2_existing_spacers_self_F, aes(x=value, y=y), color = "red") +
  scale_x_continuous(labels = human_num) +
  annotate("text", label = " II-C, self", x = 0, y = topLabelYSelf, size = fontSizeTitle, colour = "black", hjust = 0, font = globalFont)  +
  #geom_text("C1") +
  ylab("F targets (%)") +
  theme_minimal() +
  coord_cartesian(ylim = ylimit2) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

C1_self_plot_r <- ggplot(data = averageFs_self_C1_long, aes(x = Pos+binSize/2, y = value),stat="identity") + #Average data as the base data
  scale_fill_brewer(palette="Accent") +
  geom_bar(data = b_C1_self_f, stat = "identity", #Add geom_bars Ror each replicate to represent their bins
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col1, orientation = "x") +
  geom_bar(data = e_C1_self_f, stat = "identity",
           aes(x=location+binSize/2,y=target),alpha=transparency2,color=bordercolor, fill=col3, orientation = "x") +
  geom_line(data = PAMs_b185_R, stat = "density", aes(x=Minimum, y=..density..*10^7), adjust = 1/8, alpha = PAMtransparency*5, colour = "red", inherit.aes = F) + #Draw PAM Rrequency as red line
  #geom_rug(data = PAMs_b185_R, aes(x=Minimum, y=0), binwidth = 10, position = position_jitter(height = 0), sides="top", inherit.aes = F, alpha = transparency) + #This could be a spacer rug
  #geom_point(data = C2_existing_spacers_self_R, aes(x=value, y=y), color = "red") +
  geom_line(stat = 'smooth', method = 'loess', span = smoother1, size = lineThickness, alpha = 0.5, color = "black",formula = y ~ x,level = 0.999, se = FALSE) + #Draw the average as a smooth line
  ylab("R targets (%)") +
  theme_minimal() +
  #annotate("text", label = "19%", x = 1.85*10^6, y = ylimit2[2], size = annotationSize*0.7, colour = col3, alpha = 0.6)  +
  scale_y_reverse(limits=c(ylimit2[2],0)) + #Turns the reverse plot upside down
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

padding = 0.25
c1_phage_plot <- grid.arrange(C1_phage_plot_f,C1_phage_plot_r, ncol=1)
c1_self_plot <- grid.arrange(C1_self_plot_f,C1_self_plot_r, ncol=1)

c2_phage_plot <- grid.arrange(C2_phage_plot_f,C2_phage_plot_r, ncol=1)
c2_self_plot <- grid.arrange(C2_self_plot_f,C2_self_plot_r, ncol=1)

phage_plots <- grid.arrange(c1_phage_plot,c2_phage_plot,ncol=2)
self_plots <- grid.arrange(c1_self_plot, c2_self_plot, ncol=2)

uniquePlotsVertical <-grid.arrange(c1_phage_plot, c2_phage_plot, c1_self_plot, c2_self_plot, nrow = 4, ncol = 1, heights=c(padding,padding,padding,padding))
ggsave(uniquePlotsVertical, filename="fig_4B.png", height = 9.3, width = 3)

#p2 <- grid.arrange(C2_phage_plot_f,C2_self_plot_f,C2_phage_plot_r,C2_self_plot_r, ncol=2)
#### Draw plots ####
grid.newpage()

finalPositionPlot <- grid.arrange(phage_plots, self_plots)

ggsave(finalPositionPlot, filename = "fig_4B.png", height = 8.50,width = 8.90)



C1_p <- (C1_phage_plot_f / C1_phage_plot_r) + plot_layout(nrow = 2, ncol = 1)
C2_p <- (C2_phage_plot_f / C2_phage_plot_r)+ plot_layout(nrow = 2, ncol = 1)
C1_s <- (C1_self_plot_f / C1_self_plot_r)+ plot_layout(nrow = 2, ncol = 1)
C2_s <- (C2_self_plot_f / C2_self_plot_r)+ plot_layout(nrow = 2, ncol = 1)

(C1_p) / (C2_p) / (C1_s) / (C2_s) + plot_layout(nrow = 10)

