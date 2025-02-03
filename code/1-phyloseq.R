#Description -------------------------------------------------------------
#Dada2 analysis - DNA extraction comparison of 5 commercial kits
#Experiment 1 - black-capped chickadees feces in winter
#Experiment 2 - nestling blue tit feces
#Create rarefied phyloseq objects
#Sequencing: Exp1. illumina mi-seq november 2021
#            Exp 2 illumina mi-seq December 2022
#Author: Hélène Dion-Phénix
#last edition: 2023-12-04

# Libraries --------------------------------------------------------------

library(phyloseq); packageVersion("phyloseq")
library(picante); packageVersion("picante")
library(dplyr)
library(plyr)

# Import data -------------------------------------------------------------

rm(list = ls())

# ASV table
load("data/dada2/0-dada2.RData")

#metadata
load("data/0-metadata.RData")

# Phyloseq ----------------------------------------------------------------

#### create community, taxonomy, metadata files and match them by sample name
comm <- seqtab.nochim
taxo <- taxa
#### load metadata
rownames(data) <- data$sample.names
rownames(comm)[!rownames(comm) %in% rownames(data)]

# create phyloseq objects --------------------------------------------------

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),
               sample_data(data))


# Remove some taxa --------------------------------------------------------

#How many sequences are non-bacterial
(chloro = subset_taxa(ps, Order=="Chloroplast"))#5
sum(chloro@otu_table@.Data)#(250 sequences)
(euk = subset_taxa(ps, Kingdom=="Eukaryota"))#0
#(arch = subset_taxa(ps, Kingdom=="Archaea"))#1
#(mito = subset_taxa(ps, Family=='Mitochondria'))#11
sum(mito@otu_table@.Data)#2527 sequences
sum(comm)
(250+2527)/sum(comm)#0.07%

#How many sequences are not identified
(na.o = subset_taxa(ps, is.na(Order)))#911
(na.c = subset_taxa(ps, is.na(Class)))#344
(na.p = subset_taxa(ps, is.na(Phylum)))#177
sum(na.p@otu_table@.Data)
(na.k = subset_taxa(ps, is.na(Kingdom)))#2
sum(na.p@otu_table@.Data)/sum(comm)#0.6%

#Remove eucaryote, archea, Chloroplast, mitochondria and NA at phylum level
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(chloro))))
#ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(euk))))
#ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(arch))))
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(mito))))
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(na.p))))

rm(list='chloro','euk','arch','mito','na.o',"na.c","na.p","na.k")

#Extract picante/vegan format objects
comm <- otu_table(ps)
comm <- comm@.Data
taxo <- tax_table(ps)
metadata<- sample_data(ps)

all(rownames(comm) %in% rownames(data))

data <- data[rownames(comm),]

data$nseq <- rowSums(comm)
data$rich <- specnumber(comm)

data.exp <- split(data,
                  data$exp)

# Experiment 1 ------------------------------------------------------------

ps.exp1 <- subset_samples(ps,
                          exp == "1")
comm.exp1 <- otu_table(ps.exp1)
comm.exp1 <- comm.exp1@.Data
comm.exp1 <- comm.exp1[,apply(comm.exp1,2,sum)>0]

(ps.exp1 <- phyloseq(otu_table(comm.exp1, taxa_are_rows=FALSE), 
                tax_table(taxa), 
                sample_data(data)))

rarecurve(comm.exp1, step=100, label=FALSE)
rarecurve(comm.exp1, step=100, label=TRUE, xlim = c(0, 20000))

data.exp$'1'[order(data.exp$'1'$nseq),c("sample.names","type","nseq")]


# Rarefaction threshold ---------------------------------------------------

#Comparison of rarefaction threshold
par(mfrow=c(2,3))
S <- specnumber(comm.exp1)

Srare <- rarefy(comm.exp1, 1133)
p1<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "1133 (0 sample, 7 ctrls)")
abline(0, 1)

#Rarefaction at 4965 allow to keep all samples except one that had no DNA and did not amplified
Srare <- rarefy(comm.exp1, 4965)
p1<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "4965 (1 sample, 7 ctrls)")
abline(0, 1)

#Next forth threshold exclude one sample each
#Exclude PureLink-PBS-2elu-11 that have almost no band
Srare <- rarefy(comm.exp1, 9345)
p2<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "9345 (2 samples, 10 ctrls)")
abline(0, 1)

#Exclude PureLink-eau-2elu-3 that have a weak band
Srare <- rarefy(comm.exp1, 12037)
p3<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "12037 (3 samples, 12 ctrls)")
abline(0, 1)

#Exclude StoolNorgen-eau-1elu-11 that have a band
Srare <- rarefy(comm.exp1, 16047)
p4<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "16047 (4 samples, 13 ctrls)")
abline(0, 1)

#PureLink-eau-1elu-3 that have a weak to medium band
Srare <- rarefy(comm.exp1, 18637)
p5<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "18637 (5 samples, 14 ctrls)")
abline(0, 1)

rarecurve(comm.exp1, step=100,label=FALSE, xlim=c(0,30000))
abline(v = 1133)
abline(v = 4965)
abline(v = 9345)
abline(v = 12037)
abline(v = 16047)
abline(v = 18637)

#4965 could do it for maximisation of sample inclusion
#12037 seems to be a good compromise

# Experiment 2 ------------------------------------------------------------

ps.exp2 <- subset_samples(ps,
                          exp == "2")
comm.exp2 <- otu_table(ps.exp2)
comm.exp2 <- comm.exp2@.Data
comm.exp2 <- comm.exp2[,apply(comm.exp2,2,sum)>0]


(ps.exp2 <- phyloseq(otu_table(comm.exp2, taxa_are_rows=FALSE), 
                     tax_table(taxa), 
                     sample_data(data)))

data.exp$'2'[order(data.exp$'2'$nseq),c("sample.names","type","nseq")]

## Rarefaction threshold ----------------------------------------------------

rarecurve(comm.exp2, step=100, label=FALSE)
rarecurve(comm.exp2, step=100, label=TRUE, xlim = c(0, 10000))
rarecurve(comm.exp2, step=100, label=TRUE, xlim = c(0, 2000))
#Comparison of rarefaction threshold

S <- specnumber(comm.exp2) # observed number of species

#Rarefaction at 4900 allow to keep all samples except one that had no DNA and did not amplified
Srare <- rarefy(comm.exp2, 506)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)


#Rarefaction at 4900 allow to keep all samples except one that had no DNA and did not amplified
Srare <- rarefy(comm.exp2, 12225)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)

#12225 remove all ctrls and all PowerSoil samples
#12037 will work for experiment 1 and 2


# Save phyloseq object before rarefaction ---------------------------------

save.image("data/1-phyloseq_objects_notRarefied.RData")

# Rarefaction -------------------------------------------------------------

set.seed(1)
comm.r <- rrarefy(comm[which(apply(comm,1,sum)>=12037),], sample=12037)
comm.r <- comm.r[,apply(comm.r,2,sum)>0]
taxo.r <- taxo[colnames(comm.r),]
data.r <- data[rownames(comm.r),]
ps.r <- phyloseq(otu_table(comm.r, taxa_are_rows=FALSE),
                 tax_table(taxo.r),
                 sample_data(data.r))
rarecurve(comm.r, step=100,label=TRUE)


# Create phyloseq objects -------------------------------------------------

## 3 experiments - samples - Rarefied --------------------------------------

table(data.r$type)
#still 6 ctrls (3 PCR positive ctrl and 3 Norgen)

#Remove controls
ps.r <- subset_samples(ps.r, sample_data(ps.r)$type == "sample")
comm.r <-otu_table(ps.r)
comm.r <-  comm.r@.Data
taxo.r <- tax_table(ps.r)
data.r <- data.r[rownames(comm.r),]


## Experiment 1 - samples - rarefied ---------------------------------------

ps.r.exp1 <- subset_samples(ps.r, sample_data(ps.r)$exp == "1")
comm.r.exp1 <-otu_table(ps.r.exp1)
comm.r.exp1 <-  comm.r.exp1@.Data
taxo.r.exp1 <- tax_table(ps.r.exp1)
data.r.exp1 <- data.r[rownames(comm.r.exp1),]

#drop levels
data.r.exp1$kit <- droplevels(data.r.exp1$kit)
data.r.exp1$group <- droplevels(data.r.exp1$group)

## Experiment 2 - samples - rarefied ---------------------------------------

ps.r.exp2 <- subset_samples(ps.r, sample_data(ps.r)$exp == "2")
comm.r.exp2 <-otu_table(ps.r.exp2)
comm.r.exp2 <-  comm.r.exp2@.Data
taxo.r.exp2 <- tax_table(ps.r.exp2)
data.r.exp2 <- data.r[rownames(comm.r.exp2),]

#drop levels
data.r.exp2$kit <- droplevels(data.r.exp2$kit)
data.r.exp2$group <- droplevels(data.r.exp2$group)

# Clean and save environment ----------------------------------------------

rm(list = "comm.exp1", "comm.exp2", "comm.exp3", "data.exp", "metadata",
   "ps.exp1", "ps.exp2","ps.exp3","seqtab.nochim", "taxa", "track",
   "S", "Srare","p1", "p2", "p3", "p4", "p5")

save.image("data/1-phyloseq_objects.RData")

