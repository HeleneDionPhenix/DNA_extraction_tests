# DNA_extraction_tests

DNA extraction tests from bird gut microbiomes by Hélène Dion-Phénix

This repository contains R scripts for comparing different DNA extraction protocols used in the study of bacterial microbiota from the feces of black-capped chickadees and blue tits. We compared five commercial kits, two rinsing solutions prior to extraction, and two elution protocols. The code files are numbered and can be run sequentially, starting with the raw data, or any specific code can be run independently using intermediate *RData* objects. The general structure of the repository is given below. A README file is available for the data and the codes in the corresponding folders.

## Project Directory Structure

|─ DNA_extraction_tests.Rproj <br> 
|─ LICENSE <br>
|─ README.md <br>
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
|─ data <br>
│&nbsp; &nbsp; &nbsp; |─ README.md <br>
│&nbsp; &nbsp; &nbsp; |─ 0-metadata.RData <br>
│&nbsp; &nbsp; &nbsp; |─ 0-metadata_and_DNAquantification.csv <br>
│&nbsp; &nbsp; &nbsp; |─ 1-phyloseq_objects.RData <br>
│&nbsp; &nbsp; &nbsp; |─ 1-phyloseq_objects_notRarefied.RData <br>
│&nbsp; &nbsp; &nbsp; |─ 2-output_models_DNAyield.RData <br>
│&nbsp; &nbsp; &nbsp; |─ 3-output_model_diversity.RData <br>
│&nbsp; &nbsp; &nbsp; |─ 4-output_community_analyses.Data <br>
│&nbsp; &nbsp; &nbsp; |─ dada2 <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-taxa_dada2.txt <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- fastq_BCC (not included) <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- fastq_BT (not included) <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- SILVA138.1 (not included) <br>
|─ figure <br>
│&nbsp; &nbsp; &nbsp; |─ 1-DNAyield_Kit_comparison.png <br>
│&nbsp; &nbsp; &nbsp; |─ 3-Diversity_Kit_comparison.png <br>
│&nbsp; &nbsp; &nbsp; |─ 4_Barplot_Microbiota_Kit_comparison.png <br>
│&nbsp; &nbsp; &nbsp; |─ 5-PCoA_Microbiota_Kit_comparison.png <br>
│&nbsp; &nbsp; &nbsp; |─ S1_1-Rarefaction_Curves.png <br>
│&nbsp; &nbsp; &nbsp; |─ S1_2-BCC_Extraction_controls.png <br>
│&nbsp; &nbsp; &nbsp; |─ S1_3-BT_Extraction_controls.png <br>
│&nbsp; &nbsp; &nbsp; |─ S1_4-PCR_controls.png <br>
│&nbsp; &nbsp; &nbsp; |─ S3_1-PCoA_Without_StoolNorgen.png <br>
│&nbsp; &nbsp; &nbsp; |─ S4_1-PCoA_Species_comparison.png <br>
|─ table <br>
│&nbsp; &nbsp; &nbsp; |─ 2-BCC_PERMANOVA.docx <br>
│&nbsp; &nbsp; &nbsp; |─ S2_1-BCC_DNAyield.docx <br>
│&nbsp; &nbsp; &nbsp; |─ S2_1-BT_DNAyield.docx <br>
│&nbsp; &nbsp; &nbsp; |─ S2_2-BCC_Diversity.docx <br>
│&nbsp; &nbsp; &nbsp; |─ S2_2-BT-Diversity.docx <br>
│&nbsp; &nbsp; &nbsp; |─ S3_1-BCC_PERMANOVA_without_StoolNorgen.docx <br>
│&nbsp; &nbsp; &nbsp; |─ S4_1-PERMANOVA_Species_comparison.docx <br>
    
*The folder fastq_BCC, fastq_BT, and SILVA138.1 need to be added to the directory where specified above to run code 0-dada2.R*

*To load the fastq files : https://figshare.com/projects/DNA_extraction_of_bird_gut_microbiomes/220255*

*To load the SILVA Ribosomal RNA Gene Database Project, version 138.1: https://www.arb-silva.de/news/view/2024/07/11/silva-release-1382/*

