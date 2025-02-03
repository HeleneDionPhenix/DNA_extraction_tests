#Description ------------------------------------------------------------
#Chapter 1 - DNA extraction comparison of 5 commercial kits
#Experiment 1 : black-capped chickadees in winter
#Experiment 2 : nestling blue tit feces
#Formatting metadata
#Author: Hélène Dion-Phénix
#Last edition: 2023-12-01


## Import data -------------------------------------------------------------

rm(list = ls())
data <- read.csv2("data/0-Exp_1-2_metadata_and_DNAquantification.csv", 
                  sep = ";", header = T, dec = ",")


# Format variables --------------------------------------------------------

rownames(data) <- data$sample.names

data[c("species", 
     "type",
     "kit",
     "preparation",
     "n_elution",
     "group",
     "exp")] <- lapply(data[c("species", 
                               "type",
                               "kit",
                               "preparation",
                               "n_elution",
                               "group",
                               "exp")], as.factor)

# Save --------------------------------------------------------------------

# Save as RData to keep formatting
save.image("data/0-metadata.RData")
