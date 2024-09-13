# Impacts of host crop and fly sex on the Drosophila suzukii microbiome 

This repository provides supporting information and the codes for the following manuscript.

## Exploring the Impact of Host Crop and Sex on the Microbiome of Spotted Wing Drosophila: Core Communities and Co-occurrence Dynamics

Rishi Bhandari<sup>1</sup>, Adam Wong, Jana Lee,Alex Boyd<sup>*1</sup>, Laurie Agosto<sup>1</sup>, Kent Shelby <sup>1</sup>, Joseph Ringbauer Jr<sup>1</sup>, and David S. Kang Potnis<sup>1</sup>


<sup>1</sup> USDA Agricultural Research Service, Biological Control of Insects Research Laboratory, Research Park, 1503 S Providence, Columbia, MO 65203-3535, USA

<sup>2</sup> 

<sup>3</sup> 




Corresponding Author: Dave Kang, Dave.Kang@usda.gov





## Abstract



## Data Availability

Sequence data generated from this work have been deposited in the SRA (Sequencing Read Achieve) database under the BioProject accessions number PRJNA1149766. 

## Acknowledgements
Special thanks to Rebecca Schmidt for her efforts in “sample mailing coordination.”


## File tree for this project. It contains all data output from HPC to run within R R scripts are found within the R scripts directory

Description of R scripts

Preprocessing - run this script first since it generates input data for the rest of the scirpts
Diversity - general microbial diversity statistics
Network- generating and plotting networks.
.
└── SpermosphereMicrobiome2022
    ├── Bacteria
    │   ├── 16s_taxonomy.csv
    │   ├── Bacteria_spermosphere_CSS_112922.rds
    │   ├── Bacteria_spermosphere_nonnorm_112922.rds
    │   ├── otu_table_16s.csv
    │   ├── otus_16s.fasta
    │   ├── otus_16s.tre
    │   ├── otus_16s_midpoint.tre
    │   └── spermospheremetadata.csv
    ├── Differential_Abundance
    │   ├── differential_abund_121422.rds
    │   └── differential_abund_alloutput_121422.rds
    ├── Figures
    │   ├── Main Figures
    │   │   ├── PDF
    │   │   │   ├── Figure1_method.pdf
    │   │   │   ├── Figure2_imbibitiondiversity.pdf
    │   │   │   ├── Figure3_composition_top_20.pdf
    │   │   │   ├── Figure4_PCOA.pdf
    │   │   │   ├── Figure5_Diff_abundance2.pdf
    │   │   │   └── Figure6_Networks2.pdf
    │   │   └── PNG
    │   │       ├── Figure1_method.png
    │   │       ├── Figure2_imbibitiondiversity.png
    │   │       ├── Figure3_composition_top_20.png
    │   │       ├── Figure4_PCOA.png
    │   │       ├── Figure5_Diff_abundance2.png
    │   │       └── Figure6_Networks2.png
    │   └── SupplementalFigures
    │       ├── PDF
    │       │   ├── SupplementalFig1_PreliminaryLogCFU.pdf
    │       │   ├── SupplementalFig2_SequencingPerformance_120522.pdf
    │       │   ├── SupplementalFig3.pdf
    │       │   └── SupplementalFig4_composition_top_20_fungi.pdf
    │       └── PNG
    │           ├── SupplementalFig1_PreliminaryLogCFU.png
    │           ├── SupplementalFig2_SequencingPerformance_120522.png
    │           ├── SupplementalFig3.png
    │           └── SupplementalFig4_composition_top_20_fungi.png
    ├── Fungi
    │   ├── Fungi_CSSNorm_083022.rds
    │   ├── Fungi_spermosphere_unedited_083022.rds
    │   ├── METADATA.csv
    │   ├── OTU_Table.csv
    │   ├── fungal_taxonomy_DADA2_NBC.csv
    │   └── otus_R1.fasta
    ├── Networks
    │   ├── cotton_spieceasi_network.rds
    │   ├── netCompare_cotton_soil.rds
    │   ├── netCompare_cottonvssoybean.rds
    │   ├── netCompare_soybean_soil.rds
    │   ├── netConstruct_cottonvsoil.rds
    │   ├── netConstruct_cottonvssoybean.rds
    │   ├── netConstruct_soybeanvssoil.rds
    │   ├── netanalyse_cottonvsoil.rds
    │   ├── netanalyse_cottonvssoybean.rds
    │   ├── netanalyse_soybeanvssoil.rds
    │   ├── soil_spieceasi_network.rds
    │   └── soybean_spieceasi_network.rds
    ├── README.md
    ├── R_Scripts
    │   ├── Diversity.R
    │   ├── NetworkResultsPlot.R
    │   ├── Preprocessing.R
    │   ├── diff_abund.R
    │   └── transmission_analysis.R
    ├── SpermosphereMicrobiome2022.Rproj
    └── Tables
        ├── Main Tables
        │   ├── Table1_PERMANOVA_Prok_Split_120622.xlsx
        │   ├── Table2_CentralityNetworkHub.xlsx
        │   └── Table3_HubTaxa.xlsx
        └── Supplemental
            ├── SupplementalTable1_GlobalPERMANOVA_112922.xlsx
            ├── SupplementalTable2_diffabund.csv
            └── SupplementalTable3_NetworkProperties.xlsx

