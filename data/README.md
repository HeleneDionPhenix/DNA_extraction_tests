# DNA extraction tests - Data

## Metadata

The file *0-metadata_and_quantification.csv* associate each sample to the corresponding extraction protocol and give the DNA concentration (ng/uL) assess with a Qubit 3.0 high sensitivity fluorometer.

| Variable name | Description | Values |
| :------------ | ----------- | ------ |
| sample.names | Sample identifier |
| species | Bird species | [**BCC**: Black-capped chickadee, **BT**: Blue tit] |
| type | Sample type | [**sample**: feces sample, **ctrl**: PCR or extraction control] |
| kit | Extraction kit used | [**PowerSoil Pro**, **QuickDNA**, **PureLink**, **MagMAX**, **StoolNorgen**] |
| preparation | rinsing solution | [**PBS**, **water**] |
| n_elution | number of elution | [**1**, **2**] |
| group | kit - rinsing solution - number of elution | (ex: *PS_water_1*) |
| num | sample number whithin each group | [*1* to *15* and *B* (blank) for BCC, *1* to *5* and *B* for BT] |
| DNAconcentration | DNA concentration (ng/uL) assess with a Qubit 3.0 high sensitivity fluorometer | (ex: *262*) |
| sequenced | If sequenced or not | [**TRUE**/**FALSE**] |
| exp | experiment id | [**1** for BCC, **2** for BT] |


## Fastq files

The fastq files containing the 16S ARNr sequences to run these analysis are available here:
https://figshare.com/projects/DNA_extraction_of_bird_gut_microbiomes/220255

To assign the taxonomy, use the SILVA Ribosomal RNA Gene Database Project,
version 138.1:
https://www.arb-silva.de/news/view/2024/07/11/silva-release-1382/

## Intermediate data files

Alternatively, you can use the .RData files to run some specific 
part of the code without needing to download the raw data.

## Data Directory Structure

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