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
library(vegan)
require(phyloseq)
library(MicEco)

set.seed(8)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/1-phyloseq_objects.RData")

# Experiment 1 ------------------------------------------------------------

#distance matrix with hellinger distance
dist.comm.exp1<-vegdist(comm.r.exp1, method='hell')


## Multivariate homogeneity of groups dispersions -----------------------

#set order
level_order = c("MagMAX", "QuickDNA", "PowerSoil Pro", "StoolNorgen", "PureLink")
data.r.exp1$kit <- factor(data.r.exp1$kit, level = level_order)

# kit
(disper.exp1.kit<-betadisper(dist.comm.exp1, data.r.exp1$kit))
anova(disper.exp1.kit)
plot(disper.exp1.kit, label = T)

# rinsing solution
(disper.exp1.preparation<-betadisper(dist.comm.exp1, data.r.exp1$preparation))
anova(disper.exp1.preparation)
plot(disper.exp1.preparation)

# number of elution
(disper.exp1.n_elution<-betadisper(dist.comm.exp1, data.r.exp1$n_elution))
anova(disper.exp1.n_elution)
plot(disper.exp1.n_elution)

## PERMANOVA ---------------------------------------------------------------

#difference in composition between treatment
(perm.exp1<-adonis2(dist.comm.exp1 ~ kit + preparation + n_elution, 
                    data=data.r.exp1,
                    by = "term"))

#adjusted Rsquared
(perm.exp1.adj <- adonis_OmegaSq(perm.exp1))

# Experiment 2 ------------------------------------------------------------

#distance matrix with hellinger distance
dist.comm.exp2<-vegdist(comm.r.exp2, method='hell')

## Multivariate homogeneity of groups dispersions -----------------------

#set order
data.r.exp2$kit <- factor(data.r.exp2$kit, level = c("MagMAX", "PureLink"))

# kit
(disper.exp2.kit<-betadisper(dist.comm.exp2, data.r.exp2$kit))
anova(disper.exp2.kit)
plot(disper.exp2.kit, label = T)

## PERMANOVA ---------------------------------------------------------------

#difference in composition between treatment
(perm.exp2<-adonis2(dist.comm.exp2 ~ kit, 
                    data=data.r.exp2,
                    by = "term"))

#adjusted Rsquared
(perm.exp2.adj <- adonis_OmegaSq(perm.exp2))

# SUPPLEMENTARY MATERIAL --------------------------------------------------

## Experiment 1 without Norgen ------------------------------------------------

#Remove Norgen samples from phyloseq object exp1
ps.r.nonor <- subset_samples(ps.r.exp1,
                             kit != "StoolNorgen")
comm.r.nonor <- otu_table(ps.r.nonor)
comm.r.nonor <- comm.r.nonor@.Data
comm.r.nonor <- comm.r.nonor[,apply(comm.r.nonor,2,sum)>0]

#phyloseq object
(ps.r.nonor <- phyloseq(otu_table(comm.r.nonor, taxa_are_rows=FALSE), 
                        tax_table(taxo), 
                        sample_data(data)))

#metadata
data.r.nonor <- data.r[rownames(comm.r.nonor),]

#distance matrix
dist.comm.nonor <- vegdist(comm.r.nonor, method = "hellinger")


### Multivariate homogeneity of groups dispersions --------------------------

level_order = c("MagMAX", "QuickDNA", "PowerSoil Pro", "PureLink")
data.r.nonor$kit <- factor(data.r.nonor$kit, level = level_order)

# kit
(disper.nonor.kit<-betadisper(dist.comm.nonor, data.r.nonor$kit))
anova(disper.nonor.kit)
plot(disper.nonor.kit, label = T)

# rinsing solution
(disper.nonor.preparation<-betadisper(dist.comm.nonor, data.r.nonor$preparation))
anova(disper.nonor.preparation)
plot(disper.nonor.preparation)

# number of elution
(disper.nonor.n_elution<-betadisper(dist.comm.nonor, data.r.nonor$n_elution))
anova(disper.nonor.n_elution)
plot(disper.nonor.n_elution)

### PERMANOVA ---------------------------------------------------------------

#difference in composition between treatment
(perm.nonor <- adonis2(dist.comm.nonor ~ kit + preparation + n_elution, 
                       data=data.r.nonor,
                       by = "term"))

#adjusted Rsquared
(perm.nonor.adj <- adonis_OmegaSq(perm.nonor))

## Experiment 1 and 2 together ---------------------------------------------

ps.r.exp1_2 <- subset_samples(ps.r,
                             exp %in% c("1", "2"))
comm.r.exp1_2 <- otu_table(ps.r.exp1_2)
comm.r.exp1_2 <- comm.r.exp1_2@.Data
comm.r.exp1_2 <- comm.r.exp1_2[,apply(comm.r.exp1_2,2,sum)>0]

(ps.r.exp1_2 <- phyloseq(otu_table(comm.r.exp1_2, taxa_are_rows=FALSE), 
                        tax_table(taxo), 
                        sample_data(data)))

data.r.exp1_2 <- data.r[rownames(comm.r.exp1_2),]

#distance matrix with hellinger distance
dist.comm.exp1_2<-vegdist(comm.r.exp1_2, method='hell')

### Multivariate homogeneity of groups dispersions ------------------------

level_order = c("MagMAX", "QuickDNA", "PowerSoil Pro", "StoolNorgen", "PureLink")
data.r.exp1_2$kit <- factor(data.r.exp1_2$kit, level = level_order)

# kit
(disper.exp1_2.kit<-betadisper(dist.comm.exp1_2, data.r.exp1_2$kit))
anova(disper.exp1_2.kit)
plot(disper.exp1_2.kit, label = T)

# species
(disper.exp1_2.sp<-betadisper(dist.comm.exp1_2, data.r.exp1_2$species))
anova(disper.exp1_2.sp)
plot(disper.exp1_2.sp, label = T)


### PERMANOVA ---------------------------------------------------------------

#difference in composition between treatment
(perm.exp1_2<-adonis2(dist.comm.exp1_2 ~ species + kit + preparation + n_elution, 
                      data=data.r.exp1_2,
                      by = "term"))

#adjusted Rsquared
(perm.exp1_2.adj <- adonis_OmegaSq(perm.exp1_2))


# Save outputs ------------------------------------------------------------

save.image("data/3-output_community_analyses.RData")

