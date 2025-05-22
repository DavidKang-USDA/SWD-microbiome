
################################################       Differential abundance analysis          ############################################################


## Load necessary libraries
library(microeco)
library(file2meco)
library(magrittr) # for %>%
library(dplyr)    # also loads %>%

## Convert the rarefied phyloseq object to a microeco object
meco_dataset <- phyloseq2meco(bacteria.css.norm)
print(meco_dataset)

## Create a trans_abund object and plot abundance data
# Initialize trans_abund object for genus level
trans_abund_obj <- trans_abund$new(dataset = meco_dataset, taxrank = "genus")

# Plot abundance with top taxa and others in grey
trans_abund_obj$plot_bar(others_color = "grey70", facet = "Treatment", xtext_keep = FALSE, legend_text_italic = FALSE)

## Perform differential abundance analysis
# Initialize trans_diff object with LEfSe method
trans_diff_obj <- trans_diff$new(
  dataset = meco_dataset,
  method = "lefse",           # Differential abundance method
  group = "Crop",             # Grouping variable
  alpha = 0.01,               # Significance level
  lefse_subgroup = NULL       # No subgroup analysis
)

# Plot differential abundance results
# Plot bar chart with taxa sorted by LDA score (log10 scale)
trans_diff_obj$plot_diff_bar(threshold = 4) # Display taxa with LDA score above threshold
trans_diff_obj$plot_diff_bar(use_number = 1:15, width = 0.8) # Show top 15 taxa

# Save the differential abundance plot
ggsave("differentialabundanceCROP.pdf", width = 8, height = 6, units = "in", dpi = 700)

## Plot a cladogram of differential abundance
# Note: Requires ggtree package for visualization
trans_diff_obj$plot_diff_cladogram(
  use_taxa_num = 50,         # Number of taxa to include
  use_feature_num = 50,      # Number of features to include
  clade_label_level = 5       # Level of clade labels (5 corresponds to phylum level)
)


################################################    Core microbiome        ############################################################
# Load necessary libraries
library(RColorBrewer)
library(microbiome)
library(microbiomeutilities)
library(DT)
library(cooccur)
library(microViz)

# Load and preprocess the dataset (Here we are using the rerefied dataset)

# Convert to relative abundance
ps1.fly.rel <- microbiome::transform(bacteria.css.norm, "compositional")
print(ps1.fly.rel)

# Prune taxa with zero abundance
ps1.fly.rel2 <- prune_taxa(taxa_sums(ps1.fly.rel) > 0, ps1.fly.rel)
print(ps1.fly.rel2)

# Identify core taxa with detection threshold of 0.0001 and prevalence of 50%
core.taxa.standard <- core_members(ps1.fly.rel2, detection = 0.0001, prevalence = 50 / 100)
print(core.taxa.standard)

# Extract and display the taxonomy table for core OTUs
taxonomycoretaxa <- as.data.frame(tax_table(ps1.fly.rel2))
core_taxa_id <- subset(taxonomycoretaxa, rownames(taxonomycoretaxa) %in% core.taxa.standard)
DT::datatable(core_taxa_id)

# Total core abundance and diversity
core.abundance <- sample_sums(core(ps1.fly.rel2, detection = 0.0001, prevalence = 50 / 100))
DT::datatable(as.data.frame(core.abundance))

# Define parameters for heatmap plotting
prevalences <- seq(0.005, 1, 0.05)
detections <- 10^seq(log10(1e-4), log10(0.2), length = 10)

# Plot core microbiome
p.core <- plot_core(
  ps1.fly.rel2,
  plot.type = "heatmap",
  colours = rev(brewer.pal(5, "Spectral")),
  prevalences = prevalences,
  detections = detections,
  min.prevalence = 0.50
) + 
  xlab("Detection Threshold (Relative Abundance (%))")

print(p.core)

# Convert to best hit classification
ps1.stool.rel2.f <- microbiomeutilities::format_to_besthit(ps1.fly.rel2)

# Plot core microbiome with best hit classifications
p.core <- plot_core(
  ps1.stool.rel2.f,
  plot.type = "heatmap",
  colours = rev(brewer.pal(5, "Spectral")),
  prevalences = prevalences,
  detections = detections,
  min.prevalence = 0.50
) + 
  xlab("Detection Threshold (Relative Abundance (%))") +
  theme(axis.text.y = element_text(face = "italic"))

print(p.core)

# Save the core microbiome plot as a PDF
ggsave("Coremicrobiome2.pdf", width = 20, height = 6, units = "in", dpi = 700)

################################################   co-occurrence analysis       ############################################################

# Prepare data for co-occurrence analysis
ps1.core <- prune_taxa(taxa_names(ps1.fly.rel2) %in% core.taxa.standard, ps1.fly.rel2)
presenceAbsence.core <- tax_transform(ps1.core, "binary")
data.core <- otu_table(presenceAbsence.core)
data.core <- t(data.core)  # Transpose to have samples in columns and ASVs in rows

# Perform co-occurrence analysis
cooccur.core <- cooccur(
  data.core,
  type = "spp_site",
  thresh = TRUE,
  spp_names = TRUE
)
class(cooccur.core)

# Summary of co-occurrence analysis
summary(cooccur.core)

# Plot the co-occurrence analysis
p <- plot(cooccur.core, plotrand = TRUE, thresh = TRUE, spp_names = TRUE)
print(p)

# Save the co-occurrence plot as a PDF
ggsave("cooccurence.pdf", width = 7, height = 6, units = "in", dpi = 700)


################################################   Upset plot      ############################################################

# Load necessary libraries
library(VennDiagram)
library(UpSetR)
library(ggplot2)
library(MicrobiotaProcess)

# Create UpSet plot data
upsetda <- get_upset(obj = d16S, factorNames = "Treatment")

# Define plot dimensions
width_in_inches <- 12
height_in_inches <- 6

# Save UpSet plot as a PDF
pdf(file = "UPSETTreatment2.pdf", width = width_in_inches, height = height_in_inches)
upset(
  upsetda,
  sets = unique(as.vector(sample_data(d16S)$Treatment)),
  sets.bar.color = "#56B4E9",
  order.by = "freq",
  empty.intersections = "on"
)
dev.off()

# Save UpSet plot as a PNG file using ggsave
ggsave("upsetplotTreatment.pdf", width = width_in_inches, height = height_in_inches, units = "in", dpi = 700)

# Save UpSet plot data as a CSV file
write.csv(upsetda, "UpsetASVsbetweenTreatment.csv")

