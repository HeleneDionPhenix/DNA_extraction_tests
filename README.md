# DNA_extraction_tests

DNA extraction tests from bird gut microbiomes

## Project Directory Structure

|─ DNA_extraction_tests.Rproj  
|─ LICENSE  
|─ README.md  
|─ code  
│   |─ 0-dada2.R  
│   |─ 0-metadata_formatting.R  
│   |─ 1.1-phyloseq.R  
│   |─ 1.2-phyloseq_figures.R  
│   |─ 2.1-DNAyield_analyses.R  
│   |─ 2.2-DNAyield_figures.R  
│   |─ 2.3-DNAyield_tables.R  
│   |─ 3.1-diversity_analyses.R  
│   |─ 3.2-diversity_figures.R  
│   |─ 3.3-diversity_tables.R  
│   |─ 4.1-community_analyses.R  
│   |─ 4.2-community_figures.R  
│   |─ 4.3-community_tables.R  
|─ data  
│   |─ 0-metadata.RData  
│   |─ 0-metadata_and_DNAquantification.csv  
│   |─ 1-phyloseq_objects.RData  
│   |─ 1-phyloseq_objects_notRarefied.RData  
│   |─ 2-output_models_DNAyield.RData  
│   |─ 3-output_model_diversity.RData  
│   |─ 4-output_community_analyses.RData  
│   |─ README.md  
│   |─ dada2  
│       |- 0-dada2.RData  
│       |- 0-taxa_dada2.txt  
│       |- fastq_BCC (not included)  
│       |- fastq_BT (not included)  
│       |- SILVA138.1 (not included)  
|─ figure  
│   |─ 1-DNAyield_Kit_comparison.png  
│   |─ 3-Diversity_Kit_comparison.png  
│   |─ 4_Barplot_Microbiota_Kit_comparison.png  
│   |─ 5-PCoA_Microbiota_Kit_comparison.png  
│   |─ S1_1-Rarefaction_Curves.png  
│   |─ S1_2-BCC_Extraction_controls.png  
│   |─ S1_3-BT_Extraction_controls.png  
│   |─ S1_4-PCR_controls.png  
│   |─ S3_1-PCoA_Without_StoolNorgen.png  
│   └── S4_1-PCoA_Species_comparison.png  
|─ table  
    |─ 2-BCC_PERMANOVA.docx  
    |─ S2_1-BCC_DNAyield.docx  
    |─ S2_1-BT_DNAyield.docx  
    |─ S2_2-BCC_Diversity.docx  
    |─ S2_2-BT-Diversity.docx  
    |─ S3_1-BCC_PERMANOVA_without_StoolNorgen.docx  
    |─ S4_1-PERMANOVA_Species_comparison.docx  
    
*The folder fastq_BCC, fastq_BT, and SILVA138.1 need to be added to run code 0-dada2.R*
