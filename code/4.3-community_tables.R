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
library(flextable)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/3-output_community_analyses.RData")


# Tables ------------------------------------------------------------

## Experiment 1 ------------------------------------------------------------

perm.tab.exp1 <- data.frame(coef = c("Kits", "Rinsing solutions", "Numbers of elution", "Residuals", "Total"),
                            Df = perm.exp1$Df,
                            SumOfSqs = perm.exp1$SumOfSqs,
                            R2 = perm.exp1.adj$parOmegaSq,
                            F = perm.exp1$F,
                            P = perm.exp1$`Pr(>F)`)

(t.exp1 <- as_flextable(perm.tab.exp1))

save_as_docx(
  "Table 2. BCC - Bacterial composition" = t.exp1, 
  path = "table/Table_2.docx")


## Whithout norgen------------------------------------------------------------

perm.tab.nonor <- data.frame(coef = c("Kits", "Rinsing solutions", "Numbers of elution", "Residuals", "Total"),
                            Df = perm.nonor$Df,
                            SumOfSqs = perm.nonor$SumOfSqs,
                            R2 = perm.nonor.adj$parOmegaSq,
                            F = perm.nonor$F,
                            P = perm.nonor$`Pr(>F)`)

(t.nonor <- as_flextable(perm.tab.nonor))

save_as_docx(
  "Table C1" = t.nonor, 
  path = "table/Table_C1.docx")


## Two species ------------------------------------------------------------

perm.tab.exp1_2 <- data.frame(coef = c("Species", "Kits", "Rinsing solutions", "Numbers of elution", "Residuals", "Total"),
                            Df = perm.exp1_2$Df,
                            SumOfSqs = perm.exp1_2$SumOfSqs,
                            R2 = perm.exp1_2.adj$parOmegaSq,
                            F = perm.exp1_2$F,
                            P = perm.exp1_2$`Pr(>F)`)

(t.exp1_2 <- as_flextable(perm.tab.exp1_2))

save_as_docx(
  "Table C2" = t.exp1_2, 
  path = "table/Table_C2.docx")
