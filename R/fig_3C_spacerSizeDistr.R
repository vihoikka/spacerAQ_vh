library("stringr")
library("grid")
library(patchwork)

spacers_clustered_C1_b <- read.csv("../pipeline/SPACERS/B+P_b_+1_C1_spacers_clustered_C1.fasta", sep = ",", header = FALSE)
spacers_clustered_C1_e <- read.csv("../pipeline/SPACERS/B+P_e_+1_C1_spacers_clustered_C1.fasta", sep = ",", header = FALSE)
spacers_clustered_C2_b <- read.csv("../pipeline/SPACERS/B+P_b_+1_C2_spacers_clustered_C2.fasta", sep = ",", header = FALSE)
spacers_clustered_C2_d <- read.csv("../pipeline/SPACERS/B+P_d_+1_C2_spacers_clustered_C2.fasta", sep = ",", header = FALSE)
spacers_clustered_C2_e <- read.csv("../pipeline/SPACERS/B+P_e_+1_C2_spacers_clustered_C2.fasta", sep = ",", header = FALSE)
spacers_clustered_C1_b$id <- 0 #placeholder
spacers_clustered_C1_e$id <- 0 #placeholder
spacers_clustered_C2_b$id <- 0 #placeholder
spacers_clustered_C2_d$id <- 0 #placeholder
spacers_clustered_C2_e$id <- 0 #placeholder


spacers_clustered_C1_b <- spacers_clustered_C1_b[ grep(">", spacers_clustered_C1_b$V1, invert = TRUE),]
spacers_clustered_C1_e <- spacers_clustered_C1_e[ grep(">", spacers_clustered_C1_e$V1, invert = TRUE),]
spacers_clustered_C2_b <- spacers_clustered_C2_b[ grep(">", spacers_clustered_C2_b$V1, invert = TRUE),]
spacers_clustered_C2_d <- spacers_clustered_C2_d[ grep(">", spacers_clustered_C2_d$V1, invert = TRUE),]
spacers_clustered_C2_e <- spacers_clustered_C2_e[ grep(">", spacers_clustered_C2_e$V1, invert = TRUE),]


spacers_clustered_C1_b$length <- str_length(spacers_clustered_C1_b$V1)
spacers_clustered_C1_e$length <- str_length(spacers_clustered_C1_e$V1)
spacers_clustered_C2_b$length <- str_length(spacers_clustered_C2_b$V1)
spacers_clustered_C2_d$length <- str_length(spacers_clustered_C2_d$V1)
spacers_clustered_C2_e$length <- str_length(spacers_clustered_C2_e$V1)

par(mfrow=c(1,5))


C1_b_spacerplot <- barplot(table(spacers_clustered_C1_b$length),
                           #xlab = "Spacer length",
                           main = "C1 b",
                           ylim = c(0,600))

C1_e_spacerplot <- barplot(table(spacers_clustered_C1_e$length),
                           #xlab = "Spacer length",
                           main = "C1 e",
                           ylim = c(0,600))

#par(mfrow=c(1,3))


C2_b_spacerplot <- barplot(table(spacers_clustered_C2_b$length),
                           #xlab = "Spacer length",
                           main = "C2 b",
                           ylim = c(0,600))

C2_d_spacerplot <- barplot(table(spacers_clustered_C2_d$length),
                           #xlab = "Spacer length",
                           main = "C2 d",
                           ylim = c(0,600))

C2_e_spacerplot <- barplot(table(spacers_clustered_C2_e$length),
                           #xlab = "Spacer length",
                           main = "C2 e",
                           ylim = c(0,600))

theme_update(plot.title = element_text(hjust = 0.5))
C1_b_spacerplot <- ggplot(data = spacers_clustered_C1_b, aes(x=length)) + geom_histogram(binwidth = 1,color="black", fill="white") + ggtitle("II-C b")
C1_e_spacerplot <- ggplot(data = spacers_clustered_C1_e, aes(x=length)) + geom_histogram(binwidth = 1,color="black", fill="white") + ggtitle("II-C e") + ylab("")
C2_b_spacerplot <- ggplot(data = spacers_clustered_C2_b, aes(x=length)) + geom_histogram(binwidth = 1,color="black", fill="white") + ggtitle("VI-B b")
C2_d_spacerplot <- ggplot(data = spacers_clustered_C2_d, aes(x=length)) + geom_histogram(binwidth = 1,color="black", fill="white") + ggtitle("VI-B d") + ylab("")
C2_e_spacerplot <- ggplot(data = spacers_clustered_C2_e, aes(x=length)) + geom_histogram(binwidth = 1,color="black", fill="white") + ggtitle("VI-B e")
histograms <- (C1_b_spacerplot + C1_e_spacerplot) / (C2_b_spacerplot + C2_d_spacerplot) / (C2_e_spacerplot + plot_spacer())
ggsave(histograms, filename="fig_3C", height = 5, width = 3)
