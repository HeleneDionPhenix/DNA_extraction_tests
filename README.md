# DNA_extraction_tests

DNA extraction tests from bird gut microbiomes

## Project Directory Structure

|─ DNA_extraction_tests.Rproj <br> 
|─ LICENSE <br>
|─ README.md <br>
|─ code  
│&nbsp; &nbsp; &nbsp; |─ 0-dada2.R <br>
│&nbsp; &nbsp; &nbsp; |─ 0-metadata_formatting.R <br>
│&nbsp; &nbsp; &nbsp; |─ 1.1-phyloseq.R <br>
│&nbsp; &nbsp; &nbsp; |─ 1.2-phyloseq_figures.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 2.1-DNAyield_analyses.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 2.2-DNAyield_figures.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 2.3-DNAyield_tables.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 3.1-diversity_analyses.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 3.2-diversity_figures.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 3.3-diversity_tables.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 4.1-community_analyses.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 4.2-community_figures.R&nbsp; &nbsp; &nbsp;  
│&nbsp; &nbsp; &nbsp; |─ 4.3-community_tables.R&nbsp; &nbsp; &nbsp;  
|─ data  
│&nbsp; &nbsp; &nbsp; |─ 0-metadata.R  Data  
│&nbsp; &nbsp; &nbsp; |─ 0-metadata_and_DNAquantification.csv  
│&nbsp; &nbsp; &nbsp; |─ 1-phyloseq_objects.R  Data  
│&nbsp; &nbsp; &nbsp; |─ 1-phyloseq_objects_notRarefied.R  Data  
│&nbsp; &nbsp; &nbsp; |─ 2-output_models_DNAyield.R  Data  
│&nbsp; &nbsp; &nbsp; |─ 3-output_model_diversity.R  Data  
│&nbsp; &nbsp; &nbsp; |─ 4-output_community_analyses.R  Data  
│&nbsp; &nbsp; &nbsp; |─ README.md  
│&nbsp; &nbsp; &nbsp; |─ dada2  
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |- 0-dada2.R  Data  
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |- 0-taxa_dada2.txt  
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |- fastq_BCC (not included)  
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |- fastq_BT (not included)  
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |- SILVA138.1 (not included)  
|─ figure  
│&nbsp; &nbsp; &nbsp; |─ 1-DNAyield_Kit_comparison.png  
│&nbsp; &nbsp; &nbsp; |─ 3-Diversity_Kit_comparison.png  
│&nbsp; &nbsp; &nbsp; |─ 4_Barplot_Microbiota_Kit_comparison.png  
│&nbsp; &nbsp; &nbsp; |─ 5-PCoA_Microbiota_Kit_comparison.png  
│&nbsp; &nbsp; &nbsp; |─ S1_1-Rarefaction_Curves.png  
│&nbsp; &nbsp; &nbsp; |─ S1_2-BCC_Extraction_controls.png  
│&nbsp; &nbsp; &nbsp; |─ S1_3-BT_Extraction_controls.png  
│&nbsp; &nbsp; &nbsp; |─ S1_4-PCR_controls.png  
│&nbsp; &nbsp; &nbsp; |─ S3_1-PCoA_Without_StoolNorgen.png  
│&nbsp; &nbsp; &nbsp; └── S4_1-PCoA_Species_comparison.png  
|─ table  
&nbsp; &nbsp; &nbsp;  |─ 2-BCC_PERMANOVA.docx  
&nbsp; &nbsp; &nbsp;  |─ S2_1-BCC_DNAyield.docx  
&nbsp; &nbsp; &nbsp;  |─ S2_1-BT_DNAyield.docx  
&nbsp; &nbsp; &nbsp;  |─ S2_2-BCC_Diversity.docx  
&nbsp; &nbsp; &nbsp;  |─ S2_2-BT-Diversity.docx  
&nbsp; &nbsp; &nbsp;  |─ S3_1-BCC_PERMANOVA_without_StoolNorgen.docx  
&nbsp; &nbsp; &nbsp;  |─ S4_1-PERMANOVA_Species_comparison.docx  
    
*The folder fastq_BCC, fastq_BT, and SILVA138.1 need to be added to run code 0-dada2.R  *
