#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Compare DNA yield between kits and methods
#Sequencing: Exp 1. illumina mi-seq november 2021
#            Exp 2. illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2025-01-30

# Libraries --------------------------------------------------------------

library(dplyr)
library(effects)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/0-metadata.RData")

# Experiment 1 ------------------------------------------------------------

## Format data -------------------------------------------------------------

data.exp1 <- subset(data,
                    exp == "1" &
                      type == "sample")

table(data.exp1$kit, data.exp1$preparation, data.exp1$n_elution)

hist(data.exp1$DNAconcentration, breaks = 20)
hist(log(data.exp1$DNAconcentration), breaks = 20)

#distribution: log or gamma

#reorder kits by mean
data.exp1$kit <- factor(data.exp1$kit,
                        levels = c("MagMAX", 
                                   "QuickDNA", 
                                   "PowerSoil Pro", 
                                   "StoolNorgen", 
                                   "PureLink"))

## Model -----------------------------------------------------------------------

mod.exp1 <- glm(DNAconcentration ~ kit + preparation + n_elution,
             data = data.exp1,
             family = Gamma(link = "identity"))


par(mfrow = c(2,2))
plot(mod.exp1)
dev.off()

#R2
(r2.exp1 <- 1 - 121.03/154.68)

#LRT
(LRT.exp1 = 154.68-121.03)
mod.0 <- glm(DNAconcentration ~ 1,
                data = data.exp1,
                family = Gamma(link = "identity"))
(LRT.test.exp1 <- anova(mod.0, mod.exp1, test = "Chisq"))

# Estimates
summary(mod.exp1)
# Quick view of effects
plot(allEffects(mod.exp1))


## Compute confidence intervals for each condition ------------------------

# create predicted data set
dna.pred.exp1 <- distinct(select(data.exp1,
                                  kit,
                                  preparation,
                                  n_elution))
rownames(dna.pred.exp1) <- NULL

pred=predict(mod.exp1, 
             dna.pred.exp1,
             type = "link", 
             se.fit = T)

inv_link <- family(mod.exp1)$linkinv

dna.pred.exp1=data.frame(dna.pred.exp1, 
                          DNAconcentration = inv_link(pred$fit), 
                          lowerCI=inv_link(pred$fit-(1.96*pred$se.fit)), 
                          upperCI=inv_link(pred$fit+(1.96*pred$se.fit)))

# Experiment 2 ------------------------------------------------------------

## Format data -------------------------------------------------------------

data.exp2 <- subset(data,
                    exp == "2" &
                    type == "sample")

table(data.exp2$kit)

data.exp2$kit <- droplevels(data.exp2$kit)

hist(data.exp2$DNAconcentration, breaks = 20)
hist(log(data.exp2$DNAconcentration), breaks = 20)
#distribution gamma

#reorder kits by mean
data.exp2$kit <- factor(data.exp2$kit,
                        levels = c("MagMAX", 
                                   "QuickDNA", 
                                   "PowerSoil Pro", 
                                   "StoolNorgen", 
                                   "PureLink"))

## Model -------------------------------------------------------------------

#model
mod.exp2 <- glm(DNAconcentration ~ kit,
                data = data.exp2,
                family = Gamma(link = "identity"))

#model assumptions
par(mfrow = c(2,2))
plot(mod.exp2)

dev.off()

#R2
(r2.exp2 <- 1 - 17.625/69.352)

#LRT
(LRT.exp2 = 69.352-17.625)
mod.0 <- glm(DNAconcentration ~ 1,
                data = data.exp2,
                family = Gamma(link = "identity"))
(LRT.test.exp2 <- anova(mod.0, mod.exp2, test = "Chisq"))

## Compute confidence intervals for each condition ------------------------

# create predicted data set
dna.pred.exp2 <- expand_grid(kit = levels(data.exp2$kit))

pred <- predict(mod.exp2, 
                dna.pred.exp2, 
                allow.new.levels =T, 
                type = "link",
                re.form=NA, 
                se.fit = T)

inv_link <- family(mod.exp2)$linkinv

dna.pred.exp2 <- data.frame(dna.pred.exp2, 
                            DNAconcentration = inv_link(pred$fit), 
                            lowerCI = inv_link(pred$fit-(1.96*pred$se.fit)), 
                            upperCI = inv_link(pred$fit+(1.96*pred$se.fit)))

# Clean and save environment ----------------------------------------------

rm(mod.0, pred, inv_link)
save.image("data/2-output_models_DNAyield.RData")
