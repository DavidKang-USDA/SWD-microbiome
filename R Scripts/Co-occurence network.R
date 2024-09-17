################################################################################################################################################
################################################       Microbial co-occurence network with NetCoMi         ############################################################
################################################################################################################################################


# Load necessary libraries
library(igraph)
library(phyloseq)
library(SpiecEasi)
library(NetCoMi)
library(staRank)

# Load and preprocess the phyloseq object (use non normalized reads as netcomi has its own normalization)
d16S <- readRDS("C:/Users/Rishi.Bhandari/OneDrive - USDA/Desktop/Corvallis/16S.corvallis_PS.rds")
d16S <- subset_samples(d16S, ID != "COMOF12")

# Agglomerate to genus level
corvallis_genus <- tax_glom(d16S, taxrank = "genus")

################################################   Microbial co-occurence network for SWD   ############################################################

# Network construction for all samples
net_season <- netConstruct(
  corvallis_genus,  
  filtTax = "highestVar",
  filtTaxPar = list(highestVar = 92),
  measure = "spring",
  measurePar = list(nlambda = 100, rep.num = 100, Rmethod = "approx"),
  normMethod = "none",
  taxRank = "genus",
  zeroMethod = "none",
  sparsMethod = "none",
  dissFunc = "signed",
  verbose = 2,
  seed = 123456
)

# Save the network
saveRDS(net_season, "swd_spring_network.rds")

# Analyze the network
props_season <- netAnalyze(
  net_season, 
  centrLCC = FALSE,
  avDissIgnoreInf = TRUE,
  sPathNorm = FALSE,
  clustMethod = "cluster_fast_greedy",
  hubPar = c("degree", "eigenvector"),
  hubQuant = 0.9,
  lnormFit = TRUE,
  normDeg = FALSE,
  normBetw = FALSE,
  normClose = FALSE,
  normEigen = FALSE
)

# Print a summary of the network analysis
summary(props_season)

# Plot the network
pdf(file = "commonnetworkcorvallis2.pdf", width = 6, height = 5)
plot(
  props_season,
  sameLayout = TRUE,
  layoutGroup = 1,
  rmSingles = "inboth",
  nodeSize = "mclr",
  labelScale = FALSE,
  cexNodes = 1.5,
  cexLabels = 2.5,
  cexHubLabels = 3,
  cexTitle = 3.8,
  groupNames = c("male", "female"),
  hubBorderCol = "gray40"
)
dev.off()

# Detailed network plot with custom parameters
pdf(file = "detailed_network_plot.pdf", width = 6, height = 5)
plot(
  props_season,
  sameLayout = TRUE,
  layoutGroup = 1,
  repulsion = 0.7,
  shortenLabels = "intelligent",
  labelLength = 18,
  labelPattern = c(18, "'", 3),
  labelScale = FALSE,
  nodeFilter = "none",
  nodeFilterPar = 20,
  rmSingles = "inboth",
  nodeSize = "eigen",
  nodeSizeSpread = 3,
  nodeColor = "cluster",
  posCol = "darkgreen", 
  negCol = "red",
  colorVec = rainbow(13),
  edgeTranspLow = 0,
  edgeTranspHigh = 50,
  hubBorderWidth = 1.5,
  hubBorderCol = "black",
  cexLabels = 0.5,
  cexHubLabels = 0.7
)
dev.off()



################################################   Comparison network across male and females   ############################################################


# Network construction and analysis for male vs. female groups
corvallis_male <- subset_samples(corvallis_genus, Sex == "Male")
corvallis_female <- subset_samples(corvallis_genus, Sex == "Female")

# Network construction with subsetting
n_male <- nsamples(corvallis_male)
net_season_comparison <- netConstruct(
  data = corvallis_male, 
  data2 = corvallis_female,  
  filtTax = "highestVar",
  filtTaxPar = list(highestVar = 92),
  filtSamp = "highestFreq",
  filtSampPar = list(highestFreq = n_male),
  measure = "spring",
  measurePar = list(nlambda = 100, rep.num = 100, Rmethod = "approx"),
  normMethod = "none",
  taxRank = "genus",
  zeroMethod = "none",
  sparsMethod = "none",
  dissFunc = "signed",
  verbose = 2,
  seed = 123456
)

# Save the comparison network
saveRDS(net_season_comparison, "swd_malevsfemale_spring_network.rds")

# Analyze the comparison network
props_season_comparison <- netAnalyze(
  net_season_comparison, 
  centrLCC = FALSE,
  avDissIgnoreInf = TRUE,
  sPathNorm = FALSE,
  clustMethod = "cluster_fast_greedy",
  hubPar = c("degree", "eigenvector"),
  hubQuant = 0.9,
  lnormFit = TRUE,
  normDeg = FALSE,
  normBetw = FALSE,
  normClose = FALSE,
  normEigen = FALSE
)

# Print a summary of the comparison network analysis
summary(props_season_comparison)

# Plot the comparison network
pdf(file = "male_vs_female_network.pdf", width = 6, height = 5)
plot(
  props_season_comparison,
  sameLayout = TRUE,
  layoutGroup = 1,
  rmSingles = "inboth",
  nodeSize = "mclr",
  labelScale = FALSE,
  cexNodes = 1.5,
  cexLabels = 2.5,
  cexHubLabels = 3,
  cexTitle = 3.8,
  groupNames = c("male", "female"),
  hubBorderCol = "gray40"
)
dev.off()

# Detailed comparison network plot with custom parameters
pdf(file = "detailed_male_vs_female_network.pdf", width = 6, height = 5)
plot(
  props_season_comparison,
  sameLayout = TRUE,
  layoutGroup = 1,
  repulsion = 0.7,
  shortenLabels = "intelligent",
  labelLength = 18,
  labelPattern = c(18, "'", 3),
  labelScale = FALSE,
  nodeFilter = "none",
  nodeFilterPar = 20,
  rmSingles = "inboth",
  nodeSize = "eigen",
  nodeSizeSpread = 3,
  nodeColor = "cluster",
  posCol = "darkgreen", 
  negCol = "red",
  colorVec = rainbow(13),
  edgeTranspLow = 0,
  edgeTranspHigh = 50,
  hubBorderWidth = 1.5,
  hubBorderCol = "black",
  cexLabels = 0.5,
  cexHubLabels = 0.7
)
dev.off()

# Compare networks between groups
net_env_comp <- netCompare(
  props_season_comparison, 
  permTest = TRUE,
  lnormFit = FALSE,
  jaccQuant = 0.75,
  nPerm = 1000,
  cores = 30,
  seed = 20190101,
  adjust = "none"
)

# Compare networks with multiple testing adjustment
net_env_comp_adaptbh <- netCompare(
  props_season_comparison, 
  permTest = TRUE,
  lnormFit = FALSE,
  jaccQuant = 0.75,
  nPerm = 1000,
  cores = 30,
  seed = 20190101,
  adjust = "adaptBH",
  assoPerm = net_env_comp$assoPerm
)

# Print summaries of the network comparisons
summary(net_env_comp, pAdjust = TRUE, groupNames = c("Ambient", "Ozone"), digitsPval = 6)
summary(net_env_comp_adaptbh, pAdjust = TRUE, groupNames = c("Male", "Female"), digitsPval = 6)
