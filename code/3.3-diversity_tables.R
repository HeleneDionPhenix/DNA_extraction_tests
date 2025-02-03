#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Compare sequence diversity between kits and methods - figures
#Sequencing: Exp 1. illumina mi-seq november 2021
#            Exp 2. illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2023-12-07

# Librairies --------------------------------------------------------------

library(dplyr)
library(flextable)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/4-output_model_diversity.RData")

# Black-capped chickadee ------------------------------------------------------

# order
div.pred.exp1 <- div.pred.exp1[order(div.pred.exp1$kit, 
                                     div.pred.exp1$preparation, 
                                     div.pred.exp1$n_elution),]

# format numbers
div.pred.exp1$shannon <- format(round(div.pred.exp1$shannon,
                                      digits=2), 
                                nsmall = 2)
div.pred.exp1$CI <- paste("[",
                          format(round(div.pred.exp1$lowerCI,
                                       digits=2), 
                                 nsmall = 2),
                          ", ",
                          format(round(div.pred.exp1$upperCI,
                                       digits=2), 
                                 nsmall = 2),
                          "]",
                          sep = "")

div.pred.exp1 <- select(div.pred.exp1,
                        -lowerCI, 
                        -upperCI)

#rename
names(div.pred.exp1) <- c("Kit", "Rinsing solution", "Number of elution", "Estimate", "CI")



tab.f <- flextable(div.pred.exp1) %>%
  autofit()

tab.f

save_as_docx(
  "tableS2.2 - black-capped chickadee" = tab.f,
  path = "table/S2_2-BCC_Diversity.docx")


# Blue tit ----------------------------------------------------------------

# order
div.pred.exp2 <- div.pred.exp2[order(div.pred.exp2$kit),]

# format numbers
div.pred.exp2$shannon <- format(round(div.pred.exp2$shannon,
                                               digits=2), 
                                         nsmall = 2)
div.pred.exp2$CI <- paste("[",
                          format(round(div.pred.exp2$lowerCI,
                                       digits=2), 
                                 nsmall = 2),
                          ", ",
                          format(round(div.pred.exp2$upperCI,
                                       digits=2), 
                                 nsmall = 2),
                          "]",
                          sep = "")

div.pred.exp2 <- select(div.pred.exp2,
                        -lowerCI, 
                        -upperCI)

#rename
names(div.pred.exp2) <- c("Kit", "Estimate", "CI")



tab.f <- flextable(div.pred.exp2) %>%
  autofit()

tab.f

save_as_docx(
  "tableS2.2 - blue tit" = tab.f,
  path = "table/S2_2-BT-Diversity.docx")
