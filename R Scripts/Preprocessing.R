# Install and load necessary Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install specific Bioconductor packages and set version
BiocManager::install(version = '3.14')
BiocManager::install(c("dada2", "ShortRead", "Biostrings", "DECIPHER", "phangorn", "phyloseq"))

# Load libraries
library(dada2)        # Sequence data processing
library(ggplot2)      # Plotting
library(vegan)        # Community ecology analysis
library(phyloseq)     # Phylogenetic sequencing data management
library(tidyverse)    # Data manipulation and visualization
library(DECIPHER)     # Sequence analysis and manipulation
library(phangorn)     # Phylogenetic analysis tools
library(decontam)     # Contaminant detection
library(data.table)   # Data manipulation
library(phyloseq.extended) # Extended phyloseq functionalities
library(ape)          # Phylogenetics and evolution
library(DESeq2)       # Differential expression analysis
library(plotly)       # Interactive plots
library(philr)        # Compositional data analysis
library(adespatial)   # Spatial analysis of multivariate data
library(devtools)     # Package development
library(qiime2R)      # Import QIIME artifacts
library(MicrobeR)    # Data visualization
library(microbiome)   # Microbiome data analysis
library(microbiomeSeq) # Additional microbiome data analysis
library(pander)       # Rendering R objects to Pandoc's markdown
library(ranacapa)     # Data analysis
library(grid)         # Graphics support
library(gridExtra)    # Additional graphics support
library(knitr)        # Dynamic report generation
library(png)          # PNG handling
library(ggdendro)     # Dendrograms and tree plots
library(ggpubr)       # Publication-quality figures
library(RColorBrewer) # Color palettes
library(microbiomeutilities) # Utility functions
library(lme4)         # Linear mixed-effects models
library(MASS)         # Applied statistics functions
library(car)          # Regression analysis
library(emmeans)      # Estimated marginal means

# Set working directory and define paths
setwd("~/Corvallis")
path <- "~/Corvallis/Raw_reads"

# List and sort file names
fnFs <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)

# Plot quality profiles and save as PDFs
plot.quals <- plotQualityProfile(fnFs[4])
ggsave("qualplot_Fwd_WT.pdf", plot.quals, device="pdf")
plot.quals <- plotQualityProfile(fnRs[1])
ggsave("qualplot_Rev_WT.pdf", plot.quals, device="pdf")

# Read sample FASTQ file
read.fastq(fnFs[3])

# Define file paths for filtered sequences
filtpathF <- "~/Corvallis/Raw_reads"
filtpathR <- "~/Corvallis/Raw_reads"
filtFs <- list.files(filtpathF, pattern="_1.fq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="_2.fq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

# Check if forward and reverse files match
if (!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Set seed for reproducibility
set.seed(100)

# Define paths for filtered sequences
pathfiltered <- "~/Corvallis/"
filtFs <- file.path(pathfiltered, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(pathfiltered, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim sequences
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# Learn error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plot error rates and save as PDFs
errors_F <- plotErrors(errF, nominalQ=TRUE)
ggsave("errorplot_F.pdf", errors_F, device="pdf")
errors_R <- plotErrors(errR, nominalQ=TRUE)
ggsave("errorplot_R.pdf", errors_R, device="pdf")

# Dereplicate FASTQ files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Perform DADA2 denoising
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Save sequence table
saveRDS(seqtab, file = "~/Corvallis/seqtab.rds")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track the number of reads through each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/Corvallis/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/Corvallis/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# Save final results
saveRDS(seqtab.nochim, file = "~/Corvallis/seqtab_final.rds")
saveRDS(taxa, file = "~/Corvallis/tax_final.rds")

# Generate FASTA and TSV files for ASVs
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste(">ASV", seq_along(asv_seqs), sep="_")
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern = ">", replacement = "", x = asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


# Taxonomic assignment using an external training set
# Load the sequence table and training data
seqtab <- readRDS("~/Corvallis/seqtab_final.rds")

# Convert sequences from sequence table to DNAStringSet
dna <- DNAStringSet(getSequences(seqtab))

# Load the training set for taxonomic assignment
load("~/Corvallis/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET

# Perform taxonomic assignment
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # Use all processors for faster processing

# Define the ranks of interest for taxonomic classification
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Convert the output object of class "Taxa" to a matrix for easier manipulation
taxid <- t(sapply(ids, function(x) {
  # Match ranks to taxonomy and handle unclassified taxa
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

# Set column and row names for the taxonomic matrix
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab)

# Create a display-friendly version of the taxonomy matrix
taxa2 <- taxid
taxa.print2 <- taxa2 # Removing sequence rownames for display
rownames(taxa.print2) <- NULL

# View the first few rows of the taxonomy matrix
head(taxa.print2)

# Save the taxonomic assignment results to an RDS file
saveRDS(taxa2, "~/Corvallis/16Scorvallis_taxo2.rds")

