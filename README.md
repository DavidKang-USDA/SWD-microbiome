# Impacts of host crop and fly sex on the *Drosophila suzukii* microbiome 

This repository provides supporting information and the codes for the following manuscript.

## Exploring the Impact of Host Crop and Sex on the Microbiome of Spotted Wing Drosophila: Core Communities and Co-occurrence Dynamics

Rishi Bhandari<sup>1</sup>, Adam Wong, Jana Lee,Alex Boyd<sup>*1</sup>, Laurie Agosto<sup>1</sup>, Kent Shelby <sup>1</sup>, Joseph Ringbauer Jr<sup>1</sup>, and David S. Kang <sup>1</sup>


<sup>1</sup> USDA Agricultural Research Service, Biological Control of Insects Research Laboratory, Research Park, 1503 S Providence, Columbia, MO 65203-3535, USA

<sup>2</sup> 

<sup>3</sup> 




Corresponding Author: Dave Kang, Dave.Kang@usda.gov





## Abstract



## Data Availability

Sequence data generated from this work have been deposited in the SRA (Sequencing Read Achieve) database under the BioProject accessions number PRJNA1149766. 


## Acknowledgements
Special thanks to Rebecca Schmidt for her efforts in “sample mailing coordination.”

“This research used resources provided by the SCINet project and/or the AI Center of Excellence of the USDA Agricultural Research Service, ARS project numbers 0201-88888-003-000D and 0201-88888-002-000D.”


## Description of R scripts

1. [Preprocessing](https://github.com/DavidKang-USDA/SWD-microbiome/blob/main/R%20Scripts/Preprocessing.R) - This script is for preprocessing of the raw sequencing reads using [DADA2](https://benjjneb.github.io/dada2/). The script is used for generating statistics, filtering and trimming, dereplication, clustering, chimera removal, taxonomy assignment and getting a phyloseq object for downstream analysis. The output from this step is the input for subsequent steps.

2. [Diversity](https://github.com/DavidKang-USDA/SWD-microbiome/blob/main/R%20Scripts/Diversity%20analysis.R) - This script is to understand the general microbial diversity (alpha, beta and abundance) and their associated statistics

3. [Differential abundance and Core microbiome ](https://github.com/DavidKang-USDA/SWD-microbiome/blob/main/R%20Scripts/Differential%20abundance%2C%20core%2C%20cooccurence%2C%20%20and%20upset%20plot.R)- This script helps to understand differentially abundant taxa across various treatments and core microbial communiites in SWD. It also has script for conducting co-occurence pattern across the core taxa. 

4. [Network](https://github.com/DavidKang-USDA/SWD-microbiome/blob/main/R%20Scripts/Co-occurence%20network.R)- This script is for generating and plotting microbial co-occurence network in spotted wing Drosophila and comparing netwrok between male and female flies and their associated statistics.

   ### *Note: The statistics for network comparision requires huge memory and cores. So it is recommended to run this analysis in HPC.*
