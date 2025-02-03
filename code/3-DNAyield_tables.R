#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Compare DNA yield between kits and methods - tables
#Sequencing: Exp 1. illumina mi-seq november 2021
#            Exp 2. illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2025-01-30

# Librairies --------------------------------------------------------------

library(dplyr)
library(flextable)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/2-output_models_DNAyield.RData")

# Black-capped chickadee ------------------------------------------------------

# order
dna.pred.exp1 <- dna.pred.exp1[order(dna.pred.exp1$kit, 
                                     dna.pred.exp1$preparation, 
                                     dna.pred.exp1$n_elution),]

# format numbers
dna.pred.exp1$DNAconcentration <- format(round(dna.pred.exp1$DNAconcentration,
                                               digits=0), 
                                         nsmall = 0)
dna.pred.exp1$CI <- paste("[",
                          format(round(dna.pred.exp1$lowerCI,
                                       digits=0), 
                                 nsmall = 0),
                          ", ",
                          format(round(dna.pred.exp1$upperCI,
                                       digits=0), 
                                 nsmall = 0),
                          "]",
                          sep = "")

dna.pred.exp1 <- select(dna.pred.exp1,
                        -lowerCI, 
                        -upperCI)

#rename
names(dna.pred.exp1) <- c("Kit", "Rinsing solution", "Number of elution", "Estimate", "CI")



tab.f <- flextable(dna.pred.exp1) %>%
  autofit()

tab.f

save_as_docx(
  "tableS2.1 - black-capped chickadee" = tab.f,
  path = "table/S2_1-BCC_DNAyield.docx")

# Blue tit ----------------------------------------------------------------

# order
dna.pred.exp2 <- dna.pred.exp2[order(dna.pred.exp2$kit),]

# format numbers
dna.pred.exp2$DNAconcentration <- format(round(dna.pred.exp2$DNAconcentration,
                                               digits=0), 
                                         nsmall = 0)
dna.pred.exp2$CI <- paste("[",
                          format(round(dna.pred.exp2$lowerCI,
                                       digits=0), 
                                 nsmall = 0),
                          ", ",
                          format(round(dna.pred.exp2$upperCI,
                                       digits=0), 
                                 nsmall = 0),
                          "]",
                          sep = "")

dna.pred.exp2 <- select(dna.pred.exp2,
                        -lowerCI, 
                        -upperCI)

#rename
names(dna.pred.exp2) <- c("Kit", "Estimate", "CI")



tab.f <- flextable(dna.pred.exp2) %>%
  autofit()

tab.f

save_as_docx(
  "tableS2.1 - Blue tit" = tab.f,
  path = "table/S2_1-BT_DNAyield.docx")
