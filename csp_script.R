#############################################################################################################################
# CSP Analysis 
############################################################################################################################
# This script is for CSP plots for samples collected in Kenya 

# Clear the working environment
rm(list = ls())

# load packages

library(readxl)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(stringr)

# Load data 
allele_frequency <- read_excel("major_and_minor_substitutions.xlsx")



# Clean up the data to make one column wider and have Kisumu, Kilifi and Homabay 
# as factors while having the frequencies in one columns 

allele_frequency <-  allele_frequency%>%
  pivot_longer(cols = 2:4, 
               values_to = "Frequency", 
               names_to = "Region") %>%
  mutate(Allele=as.factor(Allele)) %>%
  mutate(Region=as.factor(Region))

#  Check the structure of the data 
str(allele_frequency)
 
# Plot flipped bar graph 

frequency_plot <-  ggplot(data=allele_frequency, aes(Allele, Frequency, fill=Region))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme_classic()+
  coord_flip()+
  scale_color_brewer(palette = "Dark2")+
  xlab("Allele") +
  ylab("Frequency (%)")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face="bold"))+
  scale_x_discrete(limits=rev(allele_frequency$Allele))+
  ylim(0,100)

frequency_plot
ggsave("frequency_plot.tiff", height = 9, width = 12, dpi = 300)

###############################################################################################################################
# Calculating nucleotide diversity 
##############################################################################################################################

# load data for nucleotide diversity 

CSP_pi <- read_excel("CSP_pi.xlsx")
AMA1_pi <- read_csv("AMA1_pi.csv")
ADSL_pi <- read_csv("ADSL_pi.csv")

# Edit nucleotide diversity file for CSP

colnames(CSP_pi)[1:2] <- c("Nucleotide_diversity", "pi")

# select first two columns 
CSP_pi <- CSP_pi[, 1:2]

# Check the structure of the data 
str(CSP_pi)


# Nucleotide diversity CSP 

csp_nucleotide_diversity <- ggplot(data=CSP_pi, aes(Nucleotide_diversity, pi))+
  geom_line(stat = "identity", color="blue", linewidth= 0.9)+
  theme_classic()+
  xlab("Nucleotide position")+
  ylab("pi")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  ggtitle("pfcsp")+
  theme(plot.title = element_text(colour = "black", face = "bold.italic", size = 12, hjust = 0.5))+
  xlim(0,1250)+
  ylim(0,0.07)

  

csp_nucleotide_diversity

# Nucleotide diversity AMA-1

AMA1_nucleotide_diversity <- ggplot(data=AMA1_pi, aes(NucleotideDistance, Pi))+
  geom_line(stat = "identity", color="blue", linewidth= 0.9)+
  theme_classic()+
  xlab("Nucleotide position")+
  ylab("pi")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  ggtitle("pfama1")+
  theme(plot.title = element_text(colour = "black", face = "bold.italic", size = 12, hjust = 0.5))+
  xlim(0,2000)+
  ylim(0,0.07)


AMA1_nucleotide_diversity


# ADSL nucleotide diversity

colnames(ADSL_pi)[1:2] <- c("NucleotideDistance", "Pi")

ADSL_nucleotide_diversity <- ggplot(data=ADSL_pi, aes(NucleotideDistance, Pi))+
  geom_line(stat = "identity", color="blue", linewidth= 0.9)+
  theme_classic()+
  xlab("Nucleotide position")+
  ylab("pi")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  ggtitle("pfadsl")+
  theme(plot.title = element_text(colour = "black", face = "bold.italic", size = 12, hjust = 0.5))+
  xlim(0,1500)+
  ylim(0,0.07)

ADSL_nucleotide_diversity

# Nucleotide diversity summary table 

nucleotide_diversity <- plot_grid(AMA1_nucleotide_diversity, 
                                  ADSL_nucleotide_diversity, 
                                  csp_nucleotide_diversity, 
                                  labels = c("A", "B", "C"), nrow = 1)
nucleotide_diversity

ggsave("nucleotide_diverstiy.tiff", width = 10, height = 6, dpi = 300)

##################################################################################################################################
# Tajima D graphs 
##################################################################################################################################

# load Tajima D data for controls 

CSP_tajima_D_value <- read_excel("CSP_Tajimas D.xlsx")
AMA_tajima_D <- read_csv("AMA_TajimaD.csv")
ADSL_tajima_D <- read_csv("ADSL_TajimaD.csv")

# Edit data files 

CSP_tajima_D_value <- CSP_tajima_D_value[, 2:3]
colnames(AMA_tajima_D)[1:2] <- c("Midpoint", "D")
colnames(ADSL_tajima_D)[1:2] <- c("Midpoint", "D")

# Check the structure of the data 
str(CSP_tajima_D_value)
str(AMA_tajima_D)

# Convert to become dataframe
AMA_tajima_D <- as.data.frame(AMA_tajima_D)

# plot the graphs 
csp_tajima_D<- ggplot(data=CSP_tajima_D_value, aes(Midpoint, D))+
  geom_line(stat = "identity", color="blue", linewidth= 0.9)+
  theme_classic()+
  xlab("Nucleotide position")+
  ylab("D")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  ggtitle("pfcsp")+
  theme(plot.title = element_text(colour = "black", face = "bold.italic", size = 12, hjust = 0.5))+
  xlim(0,1250)+
  ylim(-2,2)

csp_tajima_D

# Plot tajima D values for AMA1
AMA_tajima_d_plot<- ggplot(data=AMA_tajima_D, aes(Midpoint, D))+
  geom_line(stat = "identity", color="blue", linewidth= 0.9)+
  theme_classic()+
  xlab("Nucleotide position")+
  ylab("D")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  ggtitle("pfama1")+
  theme(plot.title = element_text(colour = "black", face = "bold.italic", size = 12, hjust = 0.5))+
  xlim(0,2000)+
  ylim(-2,2.5)

AMA_tajima_d_plot

# Plot Tajima D values for ADSL
ADSL_tajima_s_plot <- ggplot(data=ADSL_tajima_D, aes(Midpoint, D))+
  geom_line(stat = "identity", color="blue", linewidth= 0.9)+
  theme_classic()+
  xlab("Nucleotide position")+
  ylab("D")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  ggtitle("pfadsl")+
  theme(plot.title = element_text(colour = "black", face = "bold.italic", size = 12, hjust = 0.5))+
  xlim(0,1500)+
  ylim(-2,2.5)

ADSL_tajima_s_plot

# Use cowplot to have all three plots in one grid

Tajima_D_plots <-  plot_grid(AMA_tajima_d_plot, ADSL_tajima_s_plot, csp_tajima_D, 
                             labels = c("A", "B", "C"), nrow = 1)

Tajima_D_plots

# Save the output file locally
ggsave("Tajima_D_plots.tiff", height = 6, width = 12, dpi = 300)



#########################################################################################################
#CSP c-terminal -tajima D and mucleotide diversity
########################################################################################################

# load packages 
library(readxl)
library(tidyverse)

# load data 
 nucleotide_diversity_homabay <- read_excel("c_terminal_scp.xlsx", sheet = 2, col_names = F)
 nucleotide_diversity_kilifi <- read_excel("c_terminal_scp.xlsx", sheet = 4, col_names = F)
 nucleotide_diversity_kisumu <- read_excel("c_terminal_scp.xlsx", sheet = 6, col_names = F)


# Assign column names
  colnames(nucleotide_diversity_homabay) <- c("range", "midpoint", "pi")
  colnames(nucleotide_diversity_kisumu) <- c("range", "midpoint", "pi")
  colnames(nucleotide_diversity_kilifi) <- c("range", "midpoint", "pi")

  # extract relevant sequences 
  
  nucleotide_diversity_homabay <-  nucleotide_diversity_homabay %>%
    select(midpoint, pi) %>%
    mutate(position=midpoint +909) %>%
    select(position, pi)

nucleotide_diversity_kilifi <- nucleotide_diversity_kilifi %>%
  select(midpoint, pi) %>%
  mutate(position=midpoint +909) %>%
  select(position, pi)

nucleotide_diversity_kisumu <- nucleotide_diversity_kisumu %>%
  select(midpoint, pi) %>%
  mutate(position=midpoint +909) %>%
  select(position, pi)

# add an additional column indicating the individual regions where the samples come from 
nucleotide_diversity_homabay$site <- "Homabay"
nucleotide_diversity_kilifi$site <- "Kilifi"
nucleotide_diversity_kisumu$site <- "Kisumu"

# Join the three tables together 
 nucleotide_diversity_csp <- rbind(nucleotide_diversity_homabay,
                                   nucleotide_diversity_kisumu,
                                   nucleotide_diversity_kilifi)
 
 
 # Check structure of the data 
str(nucleotide_diversity_csp)

# convert sites into a factor
nucleotide_diversity_csp$site <- as.factor(nucleotide_diversity_csp$site)

# Plot the nucleotide diversity statistics
c_terminal_pi <- ggplot(data=nucleotide_diversity_csp, aes(position, pi, color=site))+
   geom_line(stat = "identity", linewidth= 1.0)+

     theme_classic()+
   xlab("Nucleotide position")+
   ylab("pi")+
   theme(axis.text = element_text(colour = "black", size = 10))+
   theme(axis.title = element_text(colour = "black", face = "bold", size = 11))+
   theme(plot.title = element_text(colour = "black", face = "bold", size = 12, hjust = 0.5))+
  scale_color_brewer(palette = "Set1")+
  ylim(0, 0.12)+
  scale_x_continuous(breaks = seq(900,1150, 50), 
                     labels = seq(900,1150, 50), 
                     limits = c(900,1150))+
  geom_vline(xintercept = c(930,981), linetype="dashed", color= "black") +
  geom_vline(xintercept = c(1053,1089), linetype="dashed", color= "black") +
  annotate(geom="text", label="Th3R", x= 1070, y= 0.120, color="black", fontface="bold")+
  annotate(geom="text", label="Th2R", x= 950, y= 0.120, color="black", fontface="bold")


c_terminal_pi

ggsave("c_terminal_pi.tiff", height = 6, width = 10, dpi = 300)







