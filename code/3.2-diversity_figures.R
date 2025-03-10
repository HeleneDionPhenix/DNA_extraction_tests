#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Compare sequence diversity between kits and methods - figures
#Sequencing: Exp 1. illumina mi-seq november 2021
#            Exp 2. illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2023-12-07

# Libraries --------------------------------------------------------------

library(dplyr)
library(ggplot2)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/3-output_model_diversity.RData")

# Figures ------------------------------------------------------------------

level_order = c("MagMAX", "QuickDNA", "PowerSoil Pro", "StoolNorgen", "PureLink")

div.pred.exp1$method <- paste(div.pred.exp1$preparation,
                               div.pred.exp1$n_elution,
                               sep = "-")
div.pred.exp1$method <- as.factor(div.pred.exp1$method)

data.r.exp1$method <- paste(data.r.exp1$preparation,
                          data.r.exp1$n_elution,
                          sep = "-")
data.r.exp1$method <- as.factor(data.r.exp1$method)

p.exp1 <- ggplot(div.pred.exp1, 
                 aes(x=factor(kit, level = level_order), y =shannon, color = method)) +
  
  geom_jitter(data=data.r.exp1,
              aes(color = method),
              alpha = 0.5,
              shape = 4,
              size = 1.5,
              position = position_dodge(width = 0.5),
              show.legend = F) +
  
  geom_point(data=div.pred.exp1,
             aes(x=factor(kit, level = level_order),
                 y =shannon,
                 color = method,
                 shape = method),
             size=1.5,
             position = position_dodge(width = 0.5))+
  geom_linerange(data=div.pred.exp1, 
                 aes(x=factor(kit, level = level_order), 
                     y =shannon, 
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
  
  ylim(2.5,6)+
  theme_bw() + 
  theme(text=element_text(family="Times New Roman", size=16)) +
  xlab("") +
  scale_x_discrete(labels = c("MagMAX", "QuickDNA", "PowerSoil", "StoolNorgen", "PureLink"))+
  ylab("Bacterial shannon diversity") +
  ggtitle("Black-capped chickadee") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title=element_text(size=12),
        legend.position = c(0.35,0.8),
        #legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.title=element_text(size=10),
        legend.text=element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size=10),
        legend.background = element_rect(fill='transparent', color = NA),
        legend.box.background = element_rect(fill='transparent', color = NA))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.exp1


## Experiment 2 ------------------------------------------------------------

p.exp2 <- ggplot(div.pred.exp2, aes(x=factor(kit, level = level_order), y = shannon)) +
  geom_point(data=data.r.exp2,
             shape = 4, 
             color = "orange",
             size = 1.5,
             alpha = 0.5) +
  geom_point(data=div.pred.exp2, 
             aes(x=factor(kit, level = level_order),
                 y=shannon),
             shape = 1, 
             color = "orange",
             size=2)+
  geom_linerange(data=div.pred.exp2, 
                 aes(x=factor(kit, level = level_order),
                     y=shannon, 
                     ymin= lowerCI, 
                     ymax=upperCI),
                 linewidth = 0.8,
                 color = "orange")+
  ylim(2.5,6)+
  theme_bw() + 
  xlab("") +
  scale_x_discrete(labels = c("MagMAX", "PureLink"))+
  ylab("")+
  ggtitle("Blue tit") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position= "none")+
  theme(text=element_text(family="Times New Roman", size=16)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.margin = margin(t = 5.5,r = 5.5,b = 22,l = 5.5))
p.exp2

# Composite figure for experiments 1 and 2 -----------------------------------

p.div <- ggarrange(p.exp1, p.exp2,
                   widths = c(6,3),
                   labels = "auto")
p.div

ggsave(
  "figure/Figure_3.png",
  p.div,
  width = 3740,
  height = 2494,
  units = "px",
  dpi = 500
)
