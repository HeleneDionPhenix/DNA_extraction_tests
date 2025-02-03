#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Compare sequence diversity between kits and methods
#Sequencing: Exp 1. illumina mi-seq november 2021
#            Exp 2. illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2025-01-30

# Libraries --------------------------------------------------------------

library(dplyr)
library(effects)
library(vegan)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/1-phyloseq_objects.RData")

# Experiment 1 ------------------------------------------------------------

# calcul diversity shannon
data.r.exp1$shannon <- vegan::diversity(comm.r.exp1)

hist(data.r.exp1$shannon)

#reorder kits
data.r.exp1$kit <- factor(data.r.exp1$kit,
                        levels = c("MagMAX", "QuickDNA", "PowerSoil Pro", "StoolNorgen", "PureLink"))

## Model -------------------------------------------------------------------

mod.exp1 <- lm(shannon ~ kit + preparation + n_elution,
          data = data.r.exp1)

#model assumptions
par(mfrow = c(2,2))
plot(mod.exp1)

dev.off()

summary(mod.exp1)
plot(allEffects(mod.exp1))

## Compute confidence intervals --------------------------------------------

div.pred.exp1 <- distinct(select(data.r.exp1,
                                 kit,
                                 preparation,
                                 n_elution))
rownames(div.pred.exp1) <- NULL

pred <- predict(mod.exp1, 
                div.pred.exp1,
                se.fit = T)

div.pred.exp1 <- data.frame(div.pred.exp1, 
                            shannon = pred$fit, 
                            lowerCI=pred$fit-(1.96*pred$se.fit), 
                            upperCI=pred$fit+(1.96*pred$se.fit))

# Experiment 2 ------------------------------------------------------------

# calcul diversity shannon
data.r.exp2$shannon <- vegan::diversity(comm.r.exp2)

hist(data.r.exp2$shannon)

## Model -------------------------------------------------------------------

mod.exp2 <- lm(shannon ~ kit,
                data = data.r.exp2)

#model assumptions
par(mfrow = c(2,2))
plot(mod.exp2)

dev.off()

summary(mod.exp2)
plot(allEffects(mod.exp2))


## Compute confidence intervals --------------------------------------------

# create predicted data set - shannon
div.pred.exp2 <- expand_grid(kit = levels(data.r.exp2$kit))

pred <- predict(mod.exp2, 
                div.pred.exp2,
                se.fit = T)

div.pred.exp2 <- data.frame(div.pred.exp2, 
                            shannon = pred$fit, 
                            lowerCI=pred$fit-(1.96*pred$se.fit), 
                            upperCI=pred$fit+(1.96*pred$se.fit))

# Clean and save environment ----------------------------------------------

rm(ps, ps.r, ps.r.exp1, ps.r.exp2,
   taxo, taxo.r, taxo.r.exp1, taxo.r.exp2,
    data, data.r,
   comm, comm.r, comm.r.exp1, comm.r.exp2, pred)
save.image("data/4-output_model_diversity.RData")
