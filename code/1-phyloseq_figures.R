#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Create rarefaction curves and control composition bar plots
#Supplementary Material - Section 1
#Sequencing: Exp1. illumina mi-seq november 2021
#            Exp 2 illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2025-01-30

# Libraries --------------------------------------------------------------

library(phyloseq); packageVersion("phyloseq")
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(fantaxtic)
library(scales)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/1-phyloseq_objects_notRarefied.RData")

# Rarefaction figures -----------------------------------------------

png("figure/S1_1-Rarefaction_Curves.png", width = 600, height = 600)

par(mfrow = c(2,1))

cols <- rep("black", 77)
cols[which(rowSums(comm.exp1) < 12037)] <- "red"
cols[which(data.exp$'1'$type == "ctrl")] <- "blue"
cols[which(data.exp$'1'$num == "positive")] <- "green"

rarecurve(comm.exp1, step=100,label=FALSE, xlim=c(0,15000),
          xlab = "Number of sequences",
          ylab = "Number of species",
          col = cols)
abline(v = 12037, col = "red")
mtext(side=3, adj=0, line=2.0, 'a', cex = 1.5) 


cols <- rep("black", 20)
cols[which(rowSums(comm.exp2) < 12037)] <- "red"
cols[which(data.exp$'2'$type == "ctrl")] <- "blue"
cols[which(data.exp$'2'$num == "positive")] <- "green"

rarecurve(comm.exp2, step=100,label=FALSE, xlim=c(0,15000),
          xlab = "Number of sequences",
          ylab = "Number of species",
          col = cols)
abline(v = 12037, col = "red")
mtext(side=3, adj=0, line=2.0, 'b', cex = 1.5) 

dev.off()


# Extraction control compare to samples -----------------------------------

## Black-capped chickadee -------------------------------------------------

# Subset controls and 2 samples for each kit

sub.samples.exp1 <- c()
tmp <- filter(data.exp$'1', kit != "PCR")

for (i in unique(tmp$kit)){
  sub <- tmp[which(tmp$kit == i), "sample.names"]
  sub.samples.exp1 <- c(sub.samples.exp1, sample(sub, 2))
}
sub.samples.exp1
sub.ctrl.exp1 <- data.exp$'1'[which(data.exp$'1'$type == "ctrl" & data.exp$'1'$kit != "PCR"), "sample.names"]
sub.exp1 <- c(sub.ctrl.exp1, sub.samples.exp1)
ps.sub.exp1 <- subset_samples(ps.exp1, sample_data(ps.exp1)$sample.names %in% sub.exp1)

# Get the most abundant phyla and the most abundant families within those phyla
top_nested.exp1 <- nested_top_taxa(ps.sub.exp1,
                                   top_tax_level = "Phylum",
                                   nested_tax_level = "Family",
                                   n_top_taxa = 3, 
                                   n_nested_taxa = 3,
                                   by_proportion = FALSE)

pal.exp1 <- taxon_colours(top_nested.exp1$ps_obj, tax_level = "Phylum",
                          palette = c(Actinobacteriota = "#345995",
                                      Firmicutes = "#EAC435", 
                                      Proteobacteria = "#FB4D3D",
                                      Other = "#7C8483"))
scales::show_col(pal.exp1)

# Plot the relative abundances at two levels.
p.bp.exp1.1 <- plot_nested_bar(ps_obj = top_nested.exp1$ps_obj,
                             top_level = "Phylum",
                             nested_level = "Family",
                             palette = pal.exp1,
                             relative_abundances = TRUE) +
  coord_flip() +
  facet_wrap(~kit + type,
             nrow = 5, ncol = 2,
             scales = "free")+
  ylab('Relative abundance (%)') +
  xlab("")+
  scale_y_continuous(labels = label_percent(suffix = ""))+
  theme_classic(base_size=16) +
  ggtitle("Extraction control         Feces sample")+
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=10,l=0,b=50, unit = "pt"))

p.bp.exp1.1

# Plot the absolute abundances at two levels.
p.bp.exp1.2 <- plot_nested_bar(ps_obj = top_nested.exp1$ps_obj,
                             top_level = "Phylum",
                             nested_level = "Family",
                             palette = pal.exp1,
                             relative_abundances = FALSE) +
  coord_flip() +
  facet_wrap(~kit + type,
             nrow = 5, ncol = 2,
             scales = "free")+
  ylab('Number of sequences') +
  xlab("")+
  theme_classic(base_size=16) +
  ggtitle(".  Extraction control      Feces sample")+
  scale_y_continuous(limits = c(0, 60000),
                     labels = label_number(suffix = "K", scale = 1e-3)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=0,l=0,b=50, unit = "pt"))

p.bp.exp1.2

empty <- data.frame(x = 1:5,
                    y = 1:5)

empty$lab <- c("StoolNorgen", "QuickDNA","PureLink", "PowerSoil", "MagMAX")

p.lab.kit <- ggplot(empty, aes(x=x, y=y))+
  
  geom_text(aes(x = rep(2.5, 5),
                y = c(0.5, 1.5, 2.5, 3.5, 4.5),
                label = lab),
            size=6,
            color = "black",
            angle = 90) +
  ylim(-0.5,4.9)+
  theme_void()+
  theme(panel.background = element_rect(fill = "white", color = "white"))

p.bp.exp1.1 <- ggarrange(p.lab.kit, p.bp.exp1.1,
                   nrow = 1, ncol =2,
                   widths = c(1,20))

p.ext <- ggarrange(p.bp.exp1.1, p.bp.exp1.2,
          nrow = 1, ncol =2,
          widths = c(15,20))
p.ext

ggsave("figure/S1_2-BCC_Extraction_controls.png",
       p.ext,
       dpi = 300,
       width = 10,
       height = 10)

## Blue tits ------------------------------------------------------------------

# Subset controls and 2 samples for each kit

sub.samples.exp2 <- c()
tmp <- filter(data.exp$'2', kit != "PCR")

table(tmp$kit, tmp$num)

sub.exp2 <- filter(tmp, num %in% c("negative", "1", "2"))
ps.sub.exp2 <- subset_samples(ps.exp2, sample_data(ps.exp2)$sample.names %in% sub.exp2$sample.names)

# Get the most abundant phyla and the most abundant families within those phyla
top_nested.exp2 <- nested_top_taxa(ps.sub.exp2,
                                   top_tax_level = "Phylum",
                                   nested_tax_level = "Family",
                                   n_top_taxa = 3, 
                                   n_nested_taxa = 3,
                                   by_proportion = FALSE)

pal.exp2 <- taxon_colours(top_nested.exp2$ps_obj, tax_level = "Phylum",
                          palette = c(Actinobacteriota = "#345995",
                                      Firmicutes = "#EAC435", 
                                      Proteobacteria = "#FB4D3D",
                                      Other = "#7C8483"))
scales::show_col(pal.exp2)

# Plot the relative abundances at two levels.
p.bp.exp2.1 <- plot_nested_bar(ps_obj = top_nested.exp2$ps_obj,
                               top_level = "Phylum",
                               nested_level = "Family",
                               palette = pal.exp2,
                               relative_abundances = TRUE) +
  coord_flip() +
  facet_wrap(~kit + type,
             nrow = 5, ncol = 2,
             scales = "free")+
  ylab('Relative abundance (%)') +
  xlab("")+
  scale_y_continuous(labels = label_percent(suffix = ""))+
  theme_classic(base_size=16) +
  ggtitle("Extraction control         Feces sample")+
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=10,l=0,b=50, unit = "pt"))

p.bp.exp2.1

# Plot the absolute abundances at two levels.
p.bp.exp2.2 <- plot_nested_bar(ps_obj = top_nested.exp2$ps_obj,
                               top_level = "Phylum",
                               nested_level = "Family",
                               palette = pal.exp2,
                               relative_abundances = FALSE) +
  coord_flip() +
  facet_wrap(~kit + type,
             nrow = 5, ncol = 2,
             scales = "free")+
  ylab('Number of sequences') +
  xlab("")+
  theme_classic(base_size=16) +
  ggtitle(".  Extraction control      Feces sample")+
  scale_y_continuous(limits = c(0, 40000),
                     labels = label_number(suffix = "K", scale = 1e-3)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=0,l=0,b=50, unit = "pt"))

p.bp.exp2.2

empty <- data.frame(x = 1:3,
                    y = 1:3)

empty$lab <- c("MagMAX", "PowerSoil", "PureLink")

p.lab.kit <- ggplot(empty, aes(x=x, y=y))+
  
  geom_text(aes(x = rep(2.5, 3),
                y = c(0.5, 1.5, 2.5),
                label = lab),
            size=6,
            color = "black",
            angle = 90) +
  ylim(-0.5,2.9)+
  theme_void()+
  theme(panel.background = element_rect(fill = "white", color = "white"))

p.bp.exp2.1 <- ggarrange(p.lab.kit, p.bp.exp2.1,
                         nrow = 1, ncol =2,
                         widths = c(1,20))

p.ext <- ggarrange(p.bp.exp2.1, p.bp.exp2.2,
                   nrow = 1, ncol =2,
                   widths = c(15,25))
p.ext

ggsave("figure/S1_3-BT_Extraction_controls.png",
       p.ext,
       dpi = 300,
       width = 10,
       height = 10)


# PCR control - bar plots ---------------------------

##  PCR control - Positive and negative - absolute abundance -------------

# subset positive controls
ps.pos <- subset_samples(ps, sample_data(ps)$num == "positive")
ps.pos <- prune_taxa(taxa_sums(ps.pos) > 0, ps.pos) 

# subset negative controls
ps.neg <- subset_samples(ps, sample_data(ps)$sample.names %in% c("exp1_Helene-CTRL-PCR-neg",
                                                                 "exp2_ctrl-PCR-neg-Helene"))
ps.neg <- prune_taxa(taxa_sums(ps.neg) > 0, ps.neg) 

# merge pcr controls
ps.pcr <- merge_phyloseq(ps.pos, ps.neg)

# extract legend labels and colors
p.pcr <- plot_bar(ps.pcr, fill="Genus")
g <- ggplot_build(p.pcr)
pal <- data.frame(colours = unique(g$data[[1]]["fill"]), 
                  label = unique(g$plot$data[, g$plot$labels$fill]))
pal <- setNames(pal$fill, pal$label)

# plot number of sequences
p.pcr <- plot_bar(ps.pcr, fill="Genus") +
  scale_x_discrete(name ="", 
                   labels=c("BCC",
                            "BT",
                            "BCC",
                            "BT")) +
  facet_wrap(~num,
             scales = "free_x",
             labeller = as_labeller(c("negative" = "PCR negative control",
                                      "positive" = "PCR positive control")))+
  ylab("Number of bacterial sequences") +
  scale_fill_manual(values = pal) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = 16)

p.pcr

##  PCR control - Positive - relative abundance -------------------------------

# Mock community composition 
# https://zymoresearch.eu/products/zymobiomics-microbial-community-dna-standard
# Bacillus subtilis - 12%
# Enterococcus faecalis - 12% 
# Escherichia coli - 12%
# Lactobacillus fermentum - 12%
# Listeria monocytogenes - 12%
# Pseudomonas aeruginosa - 12%
# Salmonella enterica - 12%
# Staphylococcus aureus - 12%

# Transform in relative abundance
pos.taxo <- names(sort(taxa_sums(ps.pos), decreasing=TRUE))
pos.taxo.ps <- transform_sample_counts(ps.pos, function(OTU) OTU/sum(OTU))
pos.taxo.ps<- prune_taxa(pos.taxo, pos.taxo.ps)

# Extract legend elements
p.pos <- plot_bar(pos.taxo.ps, fill="Genus")
g <- ggplot_build(p.pos)

# Change colors to fit with first plot
pal.pos <- pal
pal.pos <- pal.pos[unique(g$plot$data[, g$plot$labels$fill])]
pal.pos[which(is.na(pal.pos))] <- "grey50"

# Plot relative abundances
p.pos <- plot_bar(pos.taxo.ps, fill="Genus") +
  scale_x_discrete(name ="", 
                   labels=c("BCC",
                            "BT")) +
  scale_fill_manual(values = pal.pos,
                    guide = "none") +
  facet_wrap(~type,
             labeller = as_labeller(c("ctrl" = "PCR positive control")))+
  ylab("Bacterial relative abundance") +
  theme_bw(base_size = 16)

p.pos

## PCR negative controls -------------------------------------------------------

# Transform in relative abundance
neg.taxo <- names(sort(taxa_sums(ps.neg), decreasing=TRUE))
neg.taxo.ps <- transform_sample_counts(ps.neg, function(OTU) OTU/sum(OTU))
neg.taxo.ps<- prune_taxa(neg.taxo, neg.taxo.ps)

# Extract legend elements
p.neg <- plot_bar(neg.taxo.ps, fill="Genus") 
g <- ggplot_build(p.neg)

# Change colors to fit with first plot
pal.neg <- pal
pal.neg <- pal.neg[unique(g$plot$data[, g$plot$labels$fill])]
pal.neg[which(is.na(pal.neg))] <- "grey50"

# Plot relative abundances
p.neg <- plot_bar(neg.taxo.ps, fill="Genus") +
  scale_x_discrete(name ="", 
                   labels=c("BCC",
                            "BT")) +
  scale_fill_manual(values = pal.neg,
                    guide = "none") +
  facet_wrap(~type,
             labeller = as_labeller(c("ctrl" = "PCR negative control")))+
  ylab("Bacterial relative abundance") +
  theme_bw(base_size = 16)

p.neg

## PCR control - composite figure --------------------------------------

p.pcr.all <- ggarrange(p.neg, p.pos, p.pcr,
                       nrow = 1, ncol = 3,
                       widths = c(2,2,6),
                       labels = "auto")

p.pcr.all

ggsave("figure/S1_4-PCR_controls.png",
       p.pcr.all,
       dpi = 300,
       width = 15,
       height = 7)
