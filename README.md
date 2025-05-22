# Impacts of host crop and fly sex on the *Drosophila suzukii* microbiome 

This repository provides supporting information and the code for the following manuscript. <br />
 <br />

## Microbiome composition and co-occurrence dynamics in wild *Drosophila suzukii* are influenced by host crop, fly sex, and sampling location 

Rishi Bhandari<sup>1</sup>, Adam Chun-Nin Wong<sup>2</sup>, Jana C. Lee<sup>3</sup>, Alex Boyd<sup>*1</sup>,  Kent Shelby <sup>1</sup>, Joseph Ringbauer Jr<sup>1</sup>, and David S. Kang <sup>1</sup>  <br />
 <br />
 <br />
<sup>1</sup> Biological Control of Insects Research Laboratory, United States Department of Agriculture, Agriculture Research Service, Columbia, MO 65203, USA

<sup>2</sup> Department of Entomology and Nematology, University of Florida, Gainesville, FL 32611

<sup>3</sup> Horticultural Crops Disease and Pest Management Research Unit, United States Department of Agriculture, Agriculture Research Service, Corvallis, OR 97330, USA <br />
 <br />
 <br />       

Short title: Impacts of host crop and fly sex on the Drosophila suzukii microbiome

Corresponding Author: David S. Kang, Dave.Kang@usda.gov



    
   ### *Note: This study is currently under review.*   <br />  
  <br />  

## Description of R scripts

1. [Preprocessing](https://github.com/DavidKang-USDA/SWD-microbiome/blob/main/R%20Scripts/Preprocessing.R)—This script preprocesses the raw sequencing reads using [DADA2](https://benjjneb.github.io/dada2/). The script generates statistics, filters and trims, dereplication, clustering, chimera removal, assigns taxonomy, and gets a phyloseq object for downstream analysis. The output from this step is the input for subsequent steps.

2. [Diversity](https://github.com/DavidKang-USDA/SWD-microbiome/blob/main/R%20Scripts/Diversity%20analysis.R) - This script is to understand the general microbial diversity (alpha, beta, and abundance) and their associated statistics

3. [Differential abundance and Core microbiome ](https://github.com/DavidKang-USDA/SWD-microbiome/blob/main/R%20Scripts/Differential%20abundance%2C%20core%2C%20cooccurence%2C%20%20and%20upset%20plot.R)- This script helps to understand differentially abundant taxa across various treatments and core microbial communities in SWD. It also has a script for conducting co-occurrence patterns across the core taxa. 

4. [Network](https://github.com/DavidKang-USDA/SWD-microbiome/blob/main/R%20Scripts/Co-occurence%20network.R)- This script is for generating and plotting microbial co-occurrence networks in spotted wing Drosophila and comparing networks between male and female flies and their associated statistics.

   #### *Note: The statistics for network comparison require huge memory and cores. So it is recommended to run this analysis in HPC.*


```

├─ README.md
├─ R Scripts
│  ├─ Co-occurence network.R
│  ├─ Differential abundance, core, cooccurence,  and upset plot.R
│  ├─ Diversity analysis.R
│  └─ Preprocessing.R
├─ Intermediate files
│  ├─ Bacterial abundance
│  │  ├─ 16S.corvallis_PS.rds
│  │  ├─ 16S_taxo2.rds
│  │  ├─ ASVs.fa
│  │  ├─ ASVs_counts.tsv
│  │  ├─ ASVs_taxonomy.tsv
│  │  ├─ corvallis_metadata.csv
│  │  └─ seqtab_final.rds
│  └─ Bacterial co-occurence network
│     ├─ swd_malevsfemale_spring_network.rds
│     └─ swd_spring_network.rds
├─ Figures
│  ├─ Main figures
│  │  ├─ Fig. 1.pdf
│  │  ├─ Fig. 2.pdf
│  │  ├─ Fig. 3.pdf
│  │  ├─ Fig. 4.pdf
│  │  └─ Fig. 5.pdf
│  └─ Supplementary figures
│     ├─ Fig. S1.pdf
│     ├─ Fig. S2.pdf
│     ├─ Fig. S3.pdf
└─ Table
   └─ Supplementary tables
      ├─ Table S1.xlsx
      ├─ Table S2.xlsx
      ├─ Table S3.xlsx
      ├─ Table S4.xlsx
      ├─ Table S5.xlsx
      └─ Table S6.xlsx
```
## Data Availability

Sequence data generated from this work have been deposited in the SRA (Sequencing Read Archive) database under the BioProject accession number [PRJNA1149766](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1149766). 


## Acknowledgements
Special thanks to Rebecca Schmidt for her efforts in “sample mailing coordination.” This research used resources provided by the SCINet project and/or the AI Center of Excellence of the USDA Agricultural Research Service, ARS project numbers 0201-88888-003-000D and 0201-88888-002-000D.

## Funding
Rishi Bhandari was supported by a postdoctoral fellowship funded by the USDA Agricultural Research Service's SCINet Program and AI Center of Excellence, ARS project numbers 0201-88888-003-000D and 0201-88888-002-000D, and administered by the Oak Ridge Institute for Science and Education (ORISE) through an interagency agreement between the U.S. Department of Energy (DOE) and the U.S. Department of Agriculture (USDA). ORISE is managed by ORAU under DOE contract number DESC0014664. All opinions expressed in this paper are the authors and do not necessarily reflect the policies and views of the USDA, DOE, or ORAU/ORISE.
