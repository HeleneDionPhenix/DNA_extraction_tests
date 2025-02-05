# DNA extraction tests - Data

The code files are numbered and can be run sequentially, starting with the raw data, or any specific code can be run independently using intermediate *RData* objects.

## Code file descriptions

### 0-dada2.R

DADA2 analysis following Kembel lab script https://github.com/kembel-lab/scripts/tree/91225181acba420f2bb01acba1aa8fa4806255a1/dada2
based on dada2 tutorial https://benjjneb.github.io/dada2/tutorial.html

### 0-metadata_formatting.R

R formatting of *.csv* file. This code generate the *0-metadata.RData* file used in different subsequent scrits.

### 1.1-phyloseq.R

Rarefaction and creation of physloseq object from ASVs tables and metadata.
For more details: https://github.com/joey711/phyloseq

#### 1.2-phyloseq_figures.R

Creation of rarefaction curves and visual comparison of bacterial composition between extraction controls and samples, and between PCR positive and negative controls.

### 2.1-DNAyield_analyses.R

Comparison of DNA yield among kits, rinsing solution, and number of elution for black-capped chickadee and blue tit feces bacterial microbiota using Gamma linear models.

#### 2.2-DNAyield_figures.R

Visual comparison of DNA yield among kits, rinsing solution, and number of elution for black-capped chickadee and blue tit feces bacterial microbiota showing model estimates with 95% confidence intervals.

#### 2.3-DNAyield_tables.R

Model estimates and 95% confidence intervals of DNA yield among kits, rinsing solution, and number of elution for black-capped chickadee and blue tit feces bacterial microbiota.

### 3.1-diversity_analyses.R

Comparison of bacterial Shannon diversity among kits, rinsing solution, and number of elution for black-capped chickadee and blue tit feces bacterial microbiota using linear models.

#### 3.2-diversity_figures.R

Visual comparison of bacterial Shannon diversity among kits, rinsing solution, and number of elution for black-capped chickadee and blue tit feces bacterial microbiota showing model estimates with 95% confidence intervals.

#### 3.3-diversity_tables.R

Model estimates and 95% confidence intervals of bacterial Shannon diversity among kits, rinsing solution, and number of elution for black-capped chickadee and blue tit feces bacterial microbiota.

### 4.1-community_analyses

Tests of multivariate homogeneity of dispersion and PERMANOVA to compare respectively the among samples homogeneity and the difference in composition among kits, rinsing solution, and number of elution for black-capped chickadee and blue tit feces bacterial microbiota.

#### 4.2-community_figures.R

Comparison of the bacterial community composition of the samples of black-capped chickadee feces and blue tit feces extracted with different commercial DNA extraction kits visualise in PCoA.

#### 4.3-community_tables.R

Proportion (adjusted) if the total variation in the composition among sampled explained by kit, rinsing solution, and number of elution.

## Code Directory Structure

|─ code <br>
│&nbsp; &nbsp; &nbsp; |─ README.md <br>
│&nbsp; &nbsp; &nbsp; |─ 0-dada2.R <br>
│&nbsp; &nbsp; &nbsp; |─ 0-metadata_formatting.R <br>
│&nbsp; &nbsp; &nbsp; |─ 1.1-phyloseq.R <br>
│&nbsp; &nbsp; &nbsp; |─ 1.2-phyloseq_figures.R <br> 
│&nbsp; &nbsp; &nbsp; |─ 2.1-DNAyield_analyses.R <br>
│&nbsp; &nbsp; &nbsp; |─ 2.2-DNAyield_figures.R <br>  
│&nbsp; &nbsp; &nbsp; |─ 2.3-DNAyield_tables.R <br>
│&nbsp; &nbsp; &nbsp; |─ 3.1-diversity_analyses.R <br>
│&nbsp; &nbsp; &nbsp; |─ 3.2-diversity_figures.R <br>
│&nbsp; &nbsp; &nbsp; |─ 3.3-diversity_tables.R <br>
│&nbsp; &nbsp; &nbsp; |─ 4.1-community_analyses.R <br> 
│&nbsp; &nbsp; &nbsp; |─ 4.2-community_figures.R <br> 
│&nbsp; &nbsp; &nbsp; |─ 4.3-community_tables.R <br>