#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Compare DNA yield between kits and methods - Figures
#Sequencing: Exp 1. illumina mi-seq november 2021
#            Exp 2. illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2024-03-15

# Libraries --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

set.seed(8)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/2-output_models_DNAyield.RData")

# Figure - DNA yield between kits -----------------------------------------

## Experiment 1 ------------------------------------------------------------

dna.pred.exp1$method <- paste(dna.pred.exp1$preparation,
                               dna.pred.exp1$n_elution,
                               sep = "-")
dna.pred.exp1$method <- as.factor(dna.pred.exp1$method)

data.exp1$method <- paste(data.exp1$preparation,
                          data.exp1$n_elution,
                          sep = "-")
data.exp1$method <- as.factor(data.exp1$method)

level_order = c("MagMAX", "QuickDNA", "PowerSoil Pro", "StoolNorgen", "PureLink")


p.exp1 <- ggplot(dna.pred.exp1, 
                 aes(x=factor(kit, level = level_order), y =DNAconcentration, color = method)) +
  
  geom_jitter(data=data.exp1,
              aes(color = method),
              alpha = 0.5,
              shape = 4,
              size = 1.2,
              position = position_dodge(width = 0.5),
              show.legend = FALSE) +
  
  geom_point(data=dna.pred.exp1,
             aes(x=factor(kit, level = level_order),
                 y=DNAconcentration,
                 color = method,
                 shape = method),
             size=1.5,
             position = position_dodge(width = 0.5))+
  geom_linerange(data=dna.pred.exp1, 
                 aes(x=factor(kit, level = level_order), 
                     y=DNAconcentration, 
                     ymin= lowerCI, 
                     ymax=upperCI,
                     color = method),
                 linewidth=0.6, 
                 position = position_dodge(width = 0.5),
                 show.legend = FALSE)+
  
  scale_color_manual(name = "Rinsing solution and elution",
                     labels = c("PBS - Once", 
                                "PBS - Twice", 
                                "Water - Once", 
                                "Water - Twice"),
                     values = c("PBS-1" = "orange",
                                "PBS-2" = "darkorange4",
                                "water-1" = "#48A9A6",
                                "water-2" = "deepskyblue4"))+
  scale_shape_manual(name = "Rinsing solution and elution",
                     labels = c("PBS - Once", 
                                "PBS - Twice", 
                                "Water - Once", 
                                "Water - Twice"),
                     values = c("PBS-1" = 1,
                                "PBS-2" = 16,
                                "water-1" = 2,
                                "water-2" = 17))+
  
  ylim(0,1500)+
  theme_bw() + 
  xlab("") +
  scale_x_discrete(labels = c("MagMAX", "QuickDNA", "PowerSoil", "StoolNorgen", "PureLink"))+
  ylab("Total DNA (ng/uL)") +
  ggtitle("Black-capped chickadee") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title=element_text(size=12),
        text=element_text(family="Times New Roman", size=16),
        legend.position = c(0.45,0.8),
        legend.key.height=unit(0.4,"cm"),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.background = element_rect(fill='transparent', color = NA),
        legend.box.background = element_rect(fill='transparent', color = NA))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.exp1

## Experiment 2 ------------------------------------------------------------

level_order = c("MagMAX", "QuickDNA", "PowerSoil Pro", "StoolNorgen", "PureLink")

p.exp2 <- ggplot(dna.pred.exp2, aes(x=factor(kit, level = level_order), y =DNAconcentration)) +
  geom_point(data=data.exp2,
              shape = 4, 
              color = "orange",
             size = 1.2,
              alpha = 0.5) +
  geom_point(data=dna.pred.exp2, 
             aes(x=factor(kit, level = level_order),
                 y=DNAconcentration),
             shape = 1, 
             color = "orange",
             size=1.5)+
  geom_linerange(data=dna.pred.exp2, 
                 aes(x=factor(kit, level = level_order),
                     y=DNAconcentration, 
                     ymin= lowerCI, 
                     ymax=upperCI),
                 linewidth = 0.6,
                 color = "orange")+
  scale_x_discrete(labels = c("MagMAX", "QuickDNA", "PowerSoil", "StoolNorgen", "PureLink"))+
  ylim(0,1500)+
  theme_bw() + 
  xlab("") +
  ylab("")+
  ggtitle("Blue tit") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position= "none")+
  theme(text=element_text(family="Times New Roman", size=16)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.exp2

# Composite figure with experiments 1 and 2 -------------------------------

p.dna <- ggarrange(p.exp1, p.exp2,
                   widths = c(6,4),
                   labels = "auto")
p.dna


# Export ------------------------------------------------------------------

# Export Figure in png
ggsave(
  "figure/1-DNAyield_Kit_comparison.png",
  p.dna,
  width = 15,
  height = 10,
  units = "cm",
  dpi = 300
)
