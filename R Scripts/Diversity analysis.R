## Let's start the analysis now

# Load the dataset
d16S <- readRDS("~/Corvallis/16S.corvallis_PS.rds")

# Inspect the data
# View the first few entries of OTU table, taxonomy table, and sample data
head(otu_table(d16S))
head(tax_table(d16S))
head(sample_data(d16S))

# Additional metadata about the dataset
phy_tree(d16S)
ntaxa(d16S)
nsamples(d16S)
sample_sums(d16S)
taxa_sums(d16S)
rank_names(d16S)
sample_variables(d16S)

################################################    Normalizing reads counts to equal sequencing depth   ############################################################

# Sequencing Depth Analysis
library(data.table)

# Create a data table of sequencing depth per sample
sdt16S <- data.table(as(sample_data(d16S), "data.frame"),
                     TotalReads = sample_sums(d16S), keep.rownames = TRUE)
sdt16S 

# View samples with the lowest total reads
d16S %>%
  sample_data() %>%
  as.data.frame() %>%
  head()

# Remove the sample with the bad sequencing quality
d16S <- subset_samples(d16S, ID != "COMOF12")
d16S

# Assess the number of reads per sample and sort them
sample_sums(d16S)
sort(sample_sums(d16S))

# Histogram of sample depth
read.depths <- sample_sums(d16S)
hist(read.depths, breaks = 20)
summary(read.depths)

# Rarefaction Curves
# Create rarefaction curves and add a line indicating the minimum sample size
ggrare(d16S, step = 10, color = "Treatment", se = FALSE) +
  geom_vline(xintercept = min(sample_sums(d16S)), color = "gray60")

# Normalize data with sub-sampling
physeq_rar <- phyloseq::rarefy_even_depth(d16S, rngseed = TRUE)

# Check the number of sequences per sample after normalization
sample_sums(physeq_rar)

# Plot rarefaction curves on normalized data
p0 <- ggrare(physeq_rar, step = 10, color = "Treatment", se = TRUE)
ggsave("rarefactioncurvenorm.pdf", p0, device = "pdf")

# Facet rarefaction curves by location
p0 + facet_wrap(~Location, ncol = 4)



################################################      Alpha Diversity Analysis          ############################################################

# Transform abundance counts to proportions
count_to_prop <- function(x) { return(100 * x / sum(x)) }
d16S_4.trans <- transform_sample_counts(physeq_rar, count_to_prop)
sample_sums(d16S_4.trans)[1:5]

# Plot alpha diversity measures (Simpson and Shannon)
p_alpha <- plot_richness(d16S_4.trans, x = "Treatment", measures = c("Simpson", "Shannon"), color = "Sex") +
  geom_point(size = 5, alpha = 0.7) +
  scale_color_manual(values = c("red", "blue", "black"))

p_alpha

# Table of alpha-diversity estimators
table_r16S <- estimate_richness(physeq_rar, split = TRUE, measures = c("Observed", "InvSimpson", "Shannon", "Chao"))

# Combine alpha diversity table with sample metadata
sdr16S <- sample_data(physeq_rar)
datar16S <- cbind(sample_data(sdr16S), table_r16S) 
write.csv(datar16S, "datar16S-alphadiv.csv")

# Histogram of alpha-diversity estimators
hist(datar16S$Shannon, main = "Shannon diversity", xlab = "", breaks = 10)
hist(datar16S$InvSimpson, main = "Simpson diversity", xlab = "", breaks = 10)
hist(datar16S$Chao, main = "Chao richness", xlab = "", breaks = 15)

# Normality test for alpha diversity measures
shapiro.test(datar16S$Shannon)
shapiro.test(datar16S$Chao1)

# Plot Chao1 and Shannon diversity
A1 <- ggplot(data = datar16S, aes(x = TRT, y = Chao1)) +
  geom_point(aes(y = Chao1, group = Sex, color = Sex), position = position_dodge(width = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_boxplot(aes(fill = Sex), width = 0.5, lwd = 0.2) +
  ggtitle("Chao1 richness") +
  scale_colour_manual(values = c("Female" = "#FF7F00", "Male" = "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
  theme_grey()

# Reorder x-axis factor for Chao1 plot
newSTorder <- c("AOR", "JOR", "CROR", "COMO")
A1$data$TRT <- as.character(A1$data$TRT)
A1$data$TRT <- factor(A1$data$TRT, levels = newSTorder)

A1

A2 <- ggplot(data = datar16S, aes(x = TRT, y = Shannon)) +
  geom_point(aes(y = Shannon, group = Sex, color = Sex), position = position_dodge(width = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_boxplot(aes(fill = Sex), width = 0.5, lwd = 0.2) +
  ggtitle("Shannon diversity") +
  scale_colour_manual(values = c("Female" = "#FF7F00", "Male" = "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom") +
  theme_grey()

# Reorder x-axis factor for Shannon plot
A2$data$TRT <- as.character(A2$data$TRT)
A2$data$TRT <- factor(A2$data$TRT, levels = newSTorder)

# Combine and save plots
ggarrange(A1, A2, labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("Alpha.pdf", width = 12, height = 4.5, units = "in", dpi = 700)

# Statistical Tests on Alpha Diversity
# Kruskal-Wallis test for Shannon diversity
datar16S$CropSex <- with(datar16S, interaction(Crop, Sex))
shannonkruskal <- kruskal.test(Shannon ~ CropSex, data = datar16S)
summary(shannonkruskal)

# Pairwise comparisons for Shannon diversity
PT <- pairwise.wilcox.test(datar16S$Shannon, datar16S$CropSex, p.adjust.method = "fdr")
PT1 <- fullPTable(PT$p.value)
write.csv(PT1, "ShannonLocationSexCrop.csv")

# Kruskal-Wallis test for Chao1 diversity
datar16S$LocationSex <- with(datar16S, interaction(Location, Sex, Crop))
Chao1kruskal <- kruskal.test(Shannon ~ LocationSex, data = datar16S)
Chao1kruskal

# Pairwise comparisons for Chao1 diversity
PT <- pairwise.wilcox.test(datar16S$Chao1, datar16S$LocationSex, p.adjust.method = "fdr")
PT1 <- fullPTable(PT$p.value)
write.csv(PT1, "ShannonLocationSexCrop.csv")

# Beta Diversity Analysis
# PCoA Ordination
d16S_5.rare_pcoa <- ordinate(
  physeq = physeq_rar,
  method = "PCoA",
  distance = "bray"
)

# Plot PCoA results
B1 <- plot_ordination(
  physeq = physeq_rar,
  ordination = d16S_5.rare_pcoa,
  color = "TRT",
  shape = "Sex",
  title = "PCoA of Bacterial Communities"
) +
  scale_color_manual(values = c("#F4A460", "#4daf4a", "#1919ff", "darkorchid3")) +
  geom_point(aes(color = TRT), alpha = 1, size = 3) +
  theme_grey()

B1
ggsave("PCoA.pdf", width = 10, height = 5, units = "in", dpi = 700)

#Let’s try an NMDS instead. For NMDS plots it’s important to set a seed since the starting positions of samples in the alogrithm is random.

#Important: if you calculate your bray-curtis distance metric “in-line” it will perform a square root transformation and Wisconsin double standardization. If you don’t want this, you can calculate your bray-curtis distance separately


set.seed(1)

# Ordinate
d16S_5.rare_nmds <- ordinate(
  physeq = physeq_rar, 
  method = "NMDS", 
  distance = "bray"
)


plot_ordination(
  physeq = physeq_rar,
  ordination = d16S_5.rare_nmds,
  color = "TRT",
  shape = "Sex",
  title = "NMDS of Cullars rotation bacterial Communities"
) + 
  scale_color_manual(values = c("#F4A460",
                                "#4daf4a", "#1919ff", "darkorchid3")
  ) +
  geom_point(aes(color = TRT), alpha = 1, size = 4)

cat("stress is:", d16S_5.rare_nmds$stress)

################################################        NMDS and Vector fitting          ############################################################

# Load necessary libraries
library(vegan)
library(ggplot2)
library(gridExtra)
library(grid)
library(ape)
library(pairwiseAdonis)

# Calculate Bray-Curtis distance and perform PCoA
dist.bc <- phyloseq::distance(physeq_rar, "bray")
pcoa <- pcoa(dist.bc)

# Visualize principal components
# Plot the first ten principal components to identify their contribution to data variability
barplot(pcoa$values$Relative_eig[1:10], names.arg = paste('PCoA,', 1:10), ylab = 'Eigenvalues')

# Plot the principal components to identify any outliers
biplot.pcoa(pcoa)

# Prepare data for PERMANOVA and beta dispersion tests
# Convert OTU table to matrix and transpose if necessary
d16S_OTUs <- as(otu_table(physeq_rar), "matrix")
if (taxa_are_rows(physeq_rar)) {
  d16S_OTUs <- t(d16S_OTUs)
}

# Convert matrix to data frame and calculate Bray-Curtis distance
OTUs_scaled <- as.data.frame(d16S_OTUs)
BC.dist <- vegdist(OTUs_scaled, distance = "bray")

# Extract sample metadata
metadata <- data.frame(sample_data(physeq_rar))

# Perform PERMANOVA
adonis_result <- vegan::adonis2(BC.dist ~ metadata$TRT * metadata$Sex, permutations = 999)
print(adonis_result)

# Perform pairwise PERMANOVA
pairwise_result <- pairwise.adonis(BC.dist, metadata$Treatment)
write.csv(pairwise_result, "pairwisePermanova.csv")

# Perform beta dispersion tests
# Test for differences in dispersion between groups for 'Treatment'
dispr_treatment <- vegan::betadisper(BC.dist, metadata$Treatment)
print(dispr_treatment)

# Plot and boxplot for dispersion results
plot(dispr_treatment, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance")
boxplot(dispr_treatment, main = "", xlab = "")

# Perform permutation test for beta dispersion
permutest(dispr_treatment)


################################################        Abundance plots          ############################################################

# Load necessary libraries
library(vegan)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(viridis)
library(paletteer)

# Define color palettes
phylum_colors <- c(
  "magenta", "#599861", "#CD9BCD", "#DA5724", "red", "orange"
)

genera_colors <- c(
  "#D1A33D", "#DA5724", "#508578", "#CD9BCD",
  "orange", "#8569D5", "#5E738F", "#8A7C64", "#599861", "#CBD588",
  "#D14285", "#652926", "#C84248", "#AD6F3B", "#673770", "#5F7FC7", "magenta", "grey", "black", "cyan"
)

# Phylum-level plots

# Agglomerate to phylum-level and rename taxa
bacteria_phylum <- phyloseq::tax_glom(physeq_rar, "family")
phyloseq::taxa_names(bacteria_phylum) <- phyloseq::tax_table(bacteria_phylum)[, "family"]

# Visualize the first few rows of the OTU table
head(phyloseq::otu_table(bacteria_phylum)[1:5, 1:5])

# Melt the data for plotting
bacteria_phylum_melted <- phyloseq::psmelt(bacteria_phylum) %>%
  filter(Abundance > 0.001) %>%  # Filter out low abundance taxa
  arrange(family)                # Sort data frame alphabetically by phylum

# Plot phylum-level abundance
ggplot(bacteria_phylum_melted, aes(x = TRT, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = 0.2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Plot relative abundance at phylum level
ggplot(bacteria_phylum_melted, aes(x = Treatment, y = Abundance, fill = family)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Bacterial Phylum Across Different Treatments") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save the phylum-level plot and data
ggsave("Phylum_plot.pdf", width = 8, height = 4.5, units = "in", dpi = 700)
write.csv(bacteria_phylum_melted, "bacterial_phylum.csv")

# Genus-level plots

# Agglomerate to genus-level and transform to relative abundance
bacteria_genus <- phyloseq::tax_glom(physeq_rar, "family")
ps_rel_abund <- phyloseq::transform_sample_counts(bacteria_genus, function(x) x / sum(x))

# Get top 15 genera by abundance
top15_genera <- names(sort(taxa_sums(bacteria_genus), decreasing = TRUE))[1:10]
genus_top15 <- prune_taxa(top15_genera, ps_rel_abund)

# Melt the data for plotting
abundance_table_genus <- phyloseq::psmelt(genus_top15)

# Plot genus-level relative abundance
ggplot(abundance_table_genus, aes(x = Treatment, y = Abundance, fill = family)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = genera_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance of Top 15 Genera \n") +
  ggtitle("Genus Composition of Bacterial Communities by Treatment") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save the genus-level plot and data
ggsave("Top15genera_plot.pdf", width = 8, height = 6, units = "in", dpi = 700)
write.csv(abundance_table_genus, "Top15genera_abundance.csv")


  