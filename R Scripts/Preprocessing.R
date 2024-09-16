################################################################################################################################################
################################################        PREPROCESSING SAMPLES WITH DADA2         ############################################################
################################################################################################################################################


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
library(dplyr)        # Offers a set of tools for data manipulation
library(tibble)       # Used for creating and working with tibbles, a modern reimagining of data frames
library(Biostrings)   # Provides tools for handling biological sequences (e.g., DNA, RNA, proteins)



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

######################################################################################################################################################################################################
######################################################################################################################################################################################################
      # Lets start with phyloseq to create a phyloseq object combining all the metadata and taxonomic information
######################################################################################################################################################################################################
######################################################################################################################################################################################################

# Import metadata from CSV file
sdata <- read.csv("C:/Users/Rishi.Bhandari/OneDrive - USDA/Desktop/Corvallis/corvallis_metadata.csv",
                   row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

# Filter metadata to include only specific crops and prepare for merging with sequencing data
sdata2 <- sdata %>% 
    rownames_to_column(var = "SampleID") %>%
    filter(Crop %in% c("Blueberry", "Elderberry", "Blackberry")) %>%
    as_tibble()

head(sdata2)

# Convert SampleID to rownames for further processing
sdata3 <- sdata2 %>%
    column_to_rownames(var = "SampleID")

head(sdata3)

# Create phyloseq object by merging metadata with sequencing data
# Uncomment and modify the following lines if needed for your analysis
# sam.new <- sample_data(sdata3)
# d16S_4 <- merge_phyloseq(d16S_3, sam.new)

# Read and prepare the phyloseq object with taxonomy and ASV data
samdata <- sample_data(sdata3)
tax16S <- readRDS("16Scorvallis_taxo2.rds")
asv16S <- readRDS("seqtab_final.rds")
phyloseq_object_all <- phyloseq(
    otu_table(asv16S, taxa_are_rows = FALSE),
    tax_table(tax16S),
    sample_data(samdata)
)

# Filter ASVs based on taxonomy to exclude certain categories
ps16S_3 <- subset_taxa(phyloseq_object_all, 
                       !is.na(domain) & 
                       !order %in% c("Chloroplast") & 
                       !family %in% c("Mitochondria") &
                       !genus %in% c("Wolbachia"))

# Rename ASVs to ensure unique identifiers
dna.16S <- Biostrings::DNAStringSet(taxa_names(ps16S_3))
names(dna.16S) <- taxa_names(ps16S_3)
ps16S_4 <- merge_phyloseq(ps16S_3, dna.16S)
taxa_names(ps16S_4) <- paste0("ASV", seq(ntaxa(ps16S_4)))

# Filter samples and taxa based on abundance criteria
ps16S_5 <- prune_samples(sample_sums(ps16S_4) >= 1000, ps16S_4)
ps16S_6 <- filter_taxa(ps16S_5, function(x) sum(x) > 0, TRUE)

# Alternatively, prune taxa based on a relative abundance threshold
total.depth <- sum(otu_table(ps16S_4))
threshold <- 5e-5 * total.depth
ps16S_5 <- prune_taxa(taxa_sums(ps16S_4) > threshold, ps16S_4)

# Count unique families and genera in the filtered data
tax_table <- tax_table(ps16S_5)
families <- tax_table[, "family"]
genera <- tax_table[, "genus"]

num_families <- length(unique(families[!is.na(families)]))
num_genera <- length(unique(genera[!is.na(genera)]))

cat("Number of unique families:", num_families, "\n")
cat("Number of unique genera:", num_genera, "\n")

# Count unique genera
unique_genera <- length(unique(tax_table[, "genus"]))
print(paste("Number of unique genera:", unique_genera))

# Sequence alignment and phylogenetic tree construction
seqs.16S.nf <- refseq(ps16S_5)
alignment.16S.nf <- AlignSeqs(DNAStringSet(seqs.16S.nf), anchor = NA)
phang.align.16S.nf <- phyDat(as(alignment.16S.nf, "matrix"), type = "DNA")

# Model test and phylogenetic tree construction
mt.16S.nf <- modelTest(phang.align.16S.nf)
dm.16S.nf <- dist.ml(phang.align.16S.nf)
treeNJ.16S.nf <- NJ(dm.16S.nf)

# Fit the tree using maximum likelihood and GTR model
fit.16S.nf <- pml(treeNJ.16S.nf, data = phang.align.16S.nf)
fit.GTR.16S.nf <- update(fit.16S.nf, k = 4, inv = 0.2)
fit.GTR.16S.nf <- optim.pml(fit.GTR.16S.nf, model = "GTR", optInv = TRUE, optGamma = TRUE,
                            rearrangement = "stochastic", control = pml.control(trace = 0))

# Save the final phyloseq object with the constructed phylogenetic tree
d16S.nf <- phyloseq(
    tax_table(ps16S_5),
    sample_data(ps16S_5),
    otu_table(ps16S_5, taxa_are_rows = FALSE),
    refseq(ps16S_5),
    phy_tree(fit.GTR.16S.nf$tree)
)
saveRDS(d16S.nf, "16S.corvallis_PS.rds")


