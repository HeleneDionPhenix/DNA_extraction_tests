#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Compare bacterial composition between kits and methods
#Sequencing: Exp 1. illumina mi-seq november 2021
#            Exp 2. illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2023-12-04

# Libraries --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggordiplots)
library(microbiome)
library(fantaxtic)
library(ggnested)
library(ggnewscale)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/4-output_community_analyses.RData")

# PCoA Figures -----------------------------------------------------------------

## Experiment 1 ------------------------------------------------------------

#% explain
disper.exp1.kit$eig[1]/sum(disper.exp1.kit$eig)*100
disper.exp1.kit$eig[2]/sum(disper.exp1.kit$eig)*100
disper.exp1.kit$eig[3]/sum(disper.exp1.kit$eig)*100

# calculate 95% confidence ellipse for observations
gg.exp1 <- gg_ordiplot(ord = disper.exp1.kit,
                       groups = data.r.exp1$kit,
                       ellipse = FALSE,
                       conf = 0.95,
                       hull = FALSE,
                       spiders = FALSE)

# calculate 95% confidence ellipse for centroids
gg.exp1.se <- gg_ordiplot(ord = disper.exp1.kit,
                       groups = data.r.exp1$kit,
                       kind = 'se',
                       ellipse = FALSE,
                       conf = 0.95,
                       hull = FALSE,
                       spiders = FALSE)

#ggplot elements
names(gg.exp1)

#extract each 95% confidence ellipse for observations
MM_ellipse <- subset(gg.exp1$df_ellipse, Group == "MagMAX")
QD_ellipse <- subset(gg.exp1$df_ellipse, Group == "QuickDNA")
PS_ellipse <- subset(gg.exp1$df_ellipse, Group == "PowerSoil Pro")
SN_ellipse <- subset(gg.exp1$df_ellipse, Group == "StoolNorgen")
PL_ellipse <- subset(gg.exp1$df_ellipse, Group == "PureLink")

#extract each 95% confidence ellipse for centroids
MM_ellipse.se <- subset(gg.exp1.se$df_ellipse, Group == "MagMAX")
QD_ellipse.se <- subset(gg.exp1.se$df_ellipse, Group == "QuickDNA")
PS_ellipse.se <- subset(gg.exp1.se$df_ellipse, Group == "PowerSoil Pro")
SN_ellipse.se <- subset(gg.exp1.se$df_ellipse, Group == "StoolNorgen")
PL_ellipse.se <- subset(gg.exp1.se$df_ellipse, Group == "PureLink")

ord.exp1 <- gg.exp1$df_ord
names(ord.exp1)[3] <- "Kit"

p.exp1 <- ggplot(ord.exp1, aes(x=x, y=y, color = Kit)) +
  
  # 95% confidence ellipses - observations
  geom_polygon(data=MM_ellipse, aes(x=x, y=y), color="transparent", 
               fill="#C8AD55", alpha=0.3, linetype="solid") +
  geom_polygon(data=QD_ellipse, aes(x=x, y=y), color="transparent", 
               fill="#CC3F0C", alpha=0.3, linetype="solid") +
  geom_polygon(data=PS_ellipse, aes(x=x, y=y), color="transparent", 
               fill="#706993", alpha=0.3, linetype="solid") +
  geom_polygon(data=SN_ellipse, aes(x=x, y=y), color="transparent", 
               fill="#107E7D", alpha=0.3, linetype="solid") +
  geom_polygon(data=PL_ellipse, aes(x=x, y=y), color="transparent", 
               fill="#044B7F", alpha=0.3, linetype="solid") +
  
  # 95% confidence ellipses - centroids
  geom_path(data=MM_ellipse.se, aes(x=x, y=y), 
            color="#C8AD55", linetype="solid") +
  geom_path(data=QD_ellipse.se, aes(x=x, y=y), 
            color="#CC3F0C", linetype="solid") +
  geom_path(data=PS_ellipse.se, aes(x=x, y=y), 
            color="#706993", linetype="solid") +
  geom_path(data=SN_ellipse.se, aes(x=x, y=y), 
            color="#107E7D", linetype="solid") +
  geom_path(data=PL_ellipse.se, aes(x=x, y=y), 
            color="#044B7F", linetype="solid") +
  
  #observations
  geom_point(aes(shape = Kit,
                 color = Kit))+
  
  scale_shape_manual(values = c("MagMAX" = 1,
                                "QuickDNA" = 4,
                                "PowerSoil Pro" = 2,
                                "StoolNorgen" = 5,
                                "PureLink" = 3),
                     labels = c("MagMAX",
                                "QuickDNA",
                                "PowerSoil",
                                "StoolNorgen",
                                "PureLink"),
                     name = "Kit")+
  
  scale_color_manual(values = c("MagMAX" = "#C8AD55",
                                "QuickDNA" ="#CC3F0C",
                                "PowerSoil Pro" ="#706993",
                                "StoolNorgen" ="#107E7D",
                                "PureLink" ="#044B7F"),
                     labels = c("MagMAX",
                                "QuickDNA",
                                "PowerSoil",
                                "StoolNorgen",
                                "PureLink"),
                     name = "Kit")+
  
  xlab("PCoA 1 (44.7%)") +
  ylab("PCoA2 (5.7%)") +
  xlim(-0.6,1)+
  ylim(-0.35,0.5)+
  theme_bw(base_size = 12) +
  ggtitle("Black-capped chickadee") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.75,0.75),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA))
p.exp1

## Experiment 2 ------------------------------------------------------------

#% explain
disper.exp2.kit$eig[1]/sum(disper.exp2.kit$eig)*100
disper.exp2.kit$eig[2]/sum(disper.exp2.kit$eig)*100
disper.exp2.kit$eig[3]/sum(disper.exp2.kit$eig)*100

# calculate 95% confidence ellipse for observations
gg.exp2 <- gg_ordiplot(ord = disper.exp2.kit,
                       groups = data.r.exp2$kit,
                       ellipse = FALSE,
                       conf = 0.95,
                       hull = FALSE,
                       spiders = FALSE)

# calculate 95% confidence ellipse for centroids
gg.exp2.se <- gg_ordiplot(ord = disper.exp2.kit,
                          groups = data.r.exp2$kit,
                          kind = "se",
                          ellipse = FALSE,
                          conf = 0.95,
                          hull = FALSE,
                          spiders = FALSE)

#ggplot elements
names(gg.exp2)

#extract each 95% confidence ellipse for observations
MM_ellipse.exp2 <- subset(gg.exp2$df_ellipse, Group == "MagMAX")
PL_ellipse.exp2 <- subset(gg.exp2$df_ellipse, Group == "PureLink")

#extract each 95% confidence ellipse for centroids
MM_ellipse.exp2.se <- subset(gg.exp2.se$df_ellipse, Group == "MagMAX")
PL_ellipse.exp2.se <- subset(gg.exp2.se$df_ellipse, Group == "PureLink")

ord.exp2 <- gg.exp2$df_ord
names(ord.exp2)[3] <- "Kit"

p.exp2 <- ggplot(ord.exp2,aes(x=x, y=y, color = Kit)) +
  
  # 95% confidence ellipses - observations
  geom_polygon(data=MM_ellipse.exp2, 
               aes(x=x, y=y), 
               color = "transparent",
               fill = "#C8AD55", 
               alpha = 0.3, 
               linetype="solid") +
  geom_polygon(data=PL_ellipse.exp2, 
               aes(x=x, y=y), 
               color = "transparent", 
               fill = "#044B7C", 
               alpha = 0.3,
               linetype="solid") +
  
  # 95% confidence ellipses - centroids
  geom_path(data=MM_ellipse.exp2.se, 
               aes(x=x, y=y), 
               color = "#C8AD55", 
               linetype="solid") +
  geom_path(data=PL_ellipse.exp2.se, 
               aes(x=x, y=y), 
               color = "#044B7C", 
               linetype="solid") +
  
  # observations
  geom_point(aes(shape = Kit,
                 color = Kit))+
  
  scale_shape_manual(values = c("MagMAX" = 1, "PureLink" = 3),
                     labels = c("MagMAX",
                                "PureLink"),
                     name = "Kit")+
  scale_color_manual(values = c("MagMAX" = "#C8AD55", "PureLink" = "#044B7C"),
                     labels = c("MagMAX",
                                "PureLink"),
                     name = "Kit")+
  
  xlab("PCoA1 (34.0%)")+
  ylab("PCoA2 (13.2%)") +
  ylim(-0.35,0.5)+
  #coord_fixed(ratio = 2.5)+
  theme_bw(base_size = 12) +
  ggtitle("Blue tit") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.75,0.8),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = NA))
p.exp2

#compsite figure with PCoA ordination for experiment 1 and 2
p.pcoa <- ggarrange(p.exp1, p.exp2,
                    labels = "auto")

ggsave(
  "figure/Figure_5.png",
  p.pcoa,
  width = 3740,
  height = 2494,
  units = "px",
  dpi = 500
)

# Bar plot figures ----------------------------------------------------------------

## Experiment 1 ------------------------------------------------------------

# Get the most abundant phyla and the most abundant families within those phyla
top_nested.exp1 <- nested_top_taxa(ps.r.exp1,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Family",
                              n_top_taxa = 3, 
                              n_nested_taxa = 3)

pal.exp1 <- taxon_colours(top_nested.exp1$ps_obj, tax_level = "Phylum",
                     palette = c(Actinobacteriota = "#345995",
                                 Firmicutes = "#EAC435", 
                                 Proteobacteria = "#FB4D3D",
                                 Other = "#7C8483"))
scales::show_col(pal.exp1)

# Plot the relative abundances at two levels.
p.bp.exp1 <- plot_nested_bar(ps_obj = top_nested.exp1$ps_obj,
                             top_level = "Phylum",
                             nested_level = "Family",
                             palette = pal.exp1) +
  coord_flip() +
  facet_wrap(~factor(kit, levels = c("MagMAX", 
                                     "QuickDNA", 
                                     "PowerSoil Pro", 
                                     "StoolNorgen", 
                                     "PureLink")),
             scales = "free_y", nrow =5,ncol=1, axes = "margins",
             labeller = as_labeller(c("MagMAX" = "MagMAX",
                                      "QuickDNA" = "QuickDNA",
                                      "PowerSoil Pro" = "PowerSoil",
                                      "StoolNorgen" = "StoolNorgen",
                                      "PureLink" = "PureLink"))) +
  ylab('ASVs Relative Abundance') +
  xlab("Sample")+
  theme_classic(base_size=14) +
  ggtitle("Black-capped chickadee") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.3, 'cm'),
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(1, 'cm'), 
    plot.margin = margin(t=10,r=0,l=0,b=50, unit = "pt"))

p.bp.exp1


## Experiment 2 ----------------------------------------------------------------

# Get the most abundant phyla and the most abundant families within those phyla
top_nested <- nested_top_taxa(ps.r.exp2,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Family",
                              n_top_taxa = 3, 
                              n_nested_taxa = 3)

pal <- taxon_colours(top_nested$ps_obj, tax_level = "Phylum",
                     palette = c(Actinobacteriota = "#345995",
                                 Chloroflexi = "#3E652E", 
                                 Proteobacteria = "#FB4D3D",
                                 Other = "#7C8483"))
scales::show_col(pal)

# Plot the relative abundances at two levels.
p.bp.exp2 <- plot_nested_bar(ps_obj = top_nested$ps_obj,
                             top_level = "Phylum",
                             nested_level = "Family",
                             palette = pal) +
  coord_flip() +
  facet_wrap(~kit, nrow = 2, ncol = 1,
             scales = "free_y",
             labeller = as_labeller(c("MagMAX" = "MagMAX",
                                      "PureLink" = "PureLink"))) +
  ylab('ASVs Relative Abundance') +
  xlab('Sample')+
  theme_classic(base_size=14) +
  ggtitle("Blue tit") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.3, 'cm'),
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(1, 'cm'), 
    plot.margin = margin(t=0,r=0,l=0,b=50, unit = "pt"))

p.bp.exp2

#Compositie figure with barplots from exp1 and 2
p.bp <- ggarrange(p.bp.exp1, p.bp.exp2,
                  ncol = 1,
                  heights = c(7,3),
                  labels = "auto")
p.bp


ggsave(
  "figure/Figure_4.png",
  p.bp,
  width = 3740,
  height = 4862,
  units = "px",
  dpi = 500
)


# Supplementary Material --------------------------------------------------

## Black-capped chickadee without StoolNorgen ------------------------------

#% explain
disper.nonor.kit$eig[1]/sum(disper.nonor.kit$eig)*100
disper.nonor.kit$eig[2]/sum(disper.nonor.kit$eig)*100
disper.nonor.kit$eig[3]/sum(disper.nonor.kit$eig)*100

# calculate 95% confidence ellipse for observations
gg.nonor <- gg_ordiplot(ord = disper.nonor.kit,
                        groups = data.r.nonor$kit,
                        ellipse = FALSE,
                        conf = 0.95,
                        hull = FALSE,
                        spiders = FALSE)

# calculate 95% confidence ellipse for centroids
gg.nonor.se <- gg_ordiplot(ord = disper.nonor.kit,
                           groups = data.r.nonor$kit,
                           kind = "se",
                           ellipse = FALSE,
                           conf = 0.95,
                           hull = FALSE,
                           spiders = FALSE)

#ggplot elements
names(gg.nonor)

#extract each 95% confidence ellipse for observations
MM_ellipse_nonor <- subset(gg.nonor$df_ellipse, Group == "MagMAX")
QD_ellipse_nonor <- subset(gg.nonor$df_ellipse, Group == "QuickDNA")
PS_ellipse_nonor <- subset(gg.nonor$df_ellipse, Group == "PowerSoil Pro")
PL_ellipse_nonor <- subset(gg.nonor$df_ellipse, Group == "PureLink")

#extract each 95% confidence ellipse for centroids
MM_ellipse_nonor.se <- subset(gg.nonor.se$df_ellipse, Group == "MagMAX")
QD_ellipse_nonor.se <- subset(gg.nonor.se$df_ellipse, Group == "QuickDNA")
PS_ellipse_nonor.se <- subset(gg.nonor.se$df_ellipse, Group == "PowerSoil Pro")
PL_ellipse_nonor.se <- subset(gg.nonor.se$df_ellipse, Group == "PureLink")

ord.nonor <- gg.nonor$df_ord
names(ord.nonor)[3] <- "Kit"

p.nonor <- ggplot(ord.nonor, aes(x=x, y=y, color = Kit)) +
  
  # 95% confidence ellipses - observations
  geom_polygon(data=MM_ellipse_nonor, aes(x=x, y=y), color="transparent", fill="#C8AD55", alpha=0.3, linetype="solid") +
  geom_polygon(data=QD_ellipse_nonor, aes(x=x, y=y), color="transparent", fill="#CC3F0C", alpha=0.3, linetype="solid") +
  geom_polygon(data=PS_ellipse_nonor, aes(x=x, y=y), color="transparent", fill="#706993", alpha=0.3, linetype="solid") +
  geom_polygon(data=PL_ellipse_nonor, aes(x=x, y=y), color="transparent", fill="#044B7F", alpha=0.3, linetype="solid") +
  
  # 95% confidence ellipses - centroids
  geom_path(data=MM_ellipse_nonor.se, aes(x=x, y=y), color="#C8AD55", alpha=0.3, linetype="solid") +
  geom_path(data=QD_ellipse_nonor.se, aes(x=x, y=y), color="#CC3F0C", alpha=0.3, linetype="solid") +
  geom_path(data=PS_ellipse_nonor.se, aes(x=x, y=y), color="#706993", alpha=0.3, linetype="solid") +
  geom_path(data=PL_ellipse_nonor.se, aes(x=x, y=y), color="#044B7F", alpha=0.3, linetype="solid") +
  
  #observations
  geom_point(aes(shape = Kit,
                 color = Kit))+
  
  scale_shape_manual(values = c("MagMAX" = 1,
                                "QuickDNA" = 4,
                                "PowerSoil Pro" = 2,
                                "PureLink" = 3),
                     labels = c("MagMAX",
                                "QuickDNA",
                                "PowerSoil",
                                "PureLink"),
                     name = "Kit")+
  
  scale_color_manual(values = c("MagMAX" = "#C8AD55",
                                "QuickDNA" ="#CC3F0C",
                                "PowerSoil Pro" ="#706993",
                                "PureLink" ="#044B7F"),
                     labels = c("MagMAX",
                                "QuickDNA",
                                "PowerSoil",
                                "PureLink"),
                     name = "Kit")+
  
  xlab("PCoA 1 (13.5%)") +
  ylab("PCoA2 (9.8%)") +
  theme_bw(base_size = 10) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.75),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA))
p.nonor

ggsave(
  "figure/Figure_C1.png",
  p.nonor,
  width = 5.92,
  height = 4.07,
  dpi = 300
)

## Black-capped chickadee and blue tit feces together ---------------------

# % explain
disper.exp1_2.sp$eig[1]/sum(disper.exp1_2.sp$eig)*100
disper.exp1_2.sp$eig[2]/sum(disper.exp1_2.sp$eig)*100
disper.exp1_2.sp$eig[3]/sum(disper.exp1_2.sp$eig)*100

data.r.exp1_2$kit <- droplevels(data.r.exp1_2$kit)

# calculate 95% confidence ellipse for observations
gg.exp1_2 <- gg_ordiplot(ord = disper.exp1_2.sp,
                         groups = data.r.exp1_2$species,
                         ellipse = FALSE,
                         conf = 0.95,
                         hull = FALSE,
                         spiders = FALSE)

# calculate 95% confidence ellipse for centroids
gg.exp1_2.se <- gg_ordiplot(ord = disper.exp1_2.sp,
                            groups = data.r.exp1_2$species,
                            kind = "se",
                            ellipse = FALSE,
                            conf = 0.95,
                            hull = FALSE,
                            spiders = FALSE)

#ggplot elements
names(gg.exp1_2)

#extract each 95% confidence ellipse for observations
BC_ellipse <- subset(gg.exp1_2$df_ellipse, Group == "BCC")
BT_ellipse <- subset(gg.exp1_2$df_ellipse, Group == "BT")

#extract each 95% confidence ellipse for centroids
BC_ellipse.se <- subset(gg.exp1_2.se$df_ellipse, Group == "BCC")
BT_ellipse.se <- subset(gg.exp1_2.se$df_ellipse, Group == "BT")

ord.exp1_2 <- gg.exp1_2$df_ord

all(rownames(ord.exp1_2) == rownames(data.r.exp1_2))
ord.exp1_2 <- cbind(ord.exp1_2, data.r.exp1_2[,"kit"])
names(ord.exp1_2)[3:4] <- c("species","kit")

p.exp1_2 <- ggplot(ord.exp1_2, aes(x=x, y=y), color = kit) +

  # 95% confidence intervals - observations
  geom_polygon(data=BC_ellipse, aes(x=x, y=y, fill = "BCC"), color = "transparent", alpha = 0.2, linetype="solid") +
  geom_polygon(data=BT_ellipse, aes(x=x, y=y, fill = "BT"), color = "transparent", alpha = 0.3, linetype="solid") +
  
  # 95% confidence intervals - centroids
  geom_path(data=BC_ellipse.se, aes(x=x, y=y), color="black",linetype="solid") +
  geom_path(data=BT_ellipse.se, aes(x=x, y=y), color="#006BA6", linetype="solid") +
  
  scale_fill_manual(name = "species", 
                    values = c("BCC" =  "black",
                               "BT" = "#006BA6"),
                    labels = c("Black-capped chickadee",
                               "Blue tit"))+
  # observations
  geom_point(data = ord.exp1_2,
             aes(shape = kit,
                 color = kit))+
  
  scale_shape_manual(values = c("MagMAX" = 1,
                                "QuickDNA" = 4,
                                "PowerSoil Pro" = 2,
                                "StoolNorgen" = 5,
                                "PureLink" = 3),
                     labels = c("MagMAX",
                                "QuickDNA",
                                "PowerSoil",
                                "StoolNorgen",
                                "PureLink"),
                     name = "Kit")+
  
  scale_color_manual(values = c("MagMAX" = "#C8AD55",
                                "QuickDNA" ="#CC3F0C",
                                "PowerSoil Pro" ="#706993",
                                "StoolNorgen" ="#107E7D",
                                "PureLink" ="#044B7F"),
                     labels = c("MagMAX",
                                "QuickDNA",
                                "PowerSoil",
                                "StoolNorgen",
                                "PureLink"),
                     name = "Kit")+
  

  xlab("PCoA 1 (32.6%)") +
  ylab("PCoA2 (25.5%)") +
  theme_bw(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA))
p.exp1_2

ggsave(
  "figure/Figure_C2.png",
  p.exp1_2,
  width = 13,
  height = 8,
  units = "cm",
  dpi = 300
)
