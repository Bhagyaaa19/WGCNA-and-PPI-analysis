# 18/09/2024 - Bhagya Wijeratne - 15219 

# script to perform WGCNA on maize dataset



#installing the packages 
install.packages("BiocManager")
BiocManager::install("GO.db",dependencies = TRUE,force = TRUE)
BiocManager::install("WGCNA",dependencies = TRUE,force = TRUE)
BiocManager::install("DESeq2",dependencies = TRUE)
BiocManager::install("GEOquery",dependencies = TRUE,force = TRUE)
install.packages(c("ggplot2", "gridExtra", "reshape2"))
BiocManager::install(c("ComplexHeatmap", "circlize"), force = TRUE)
install.packages("devtools")
devtools::install_github("kevinblighe/CorLevelPlot")
install.packages("readxl") # to read excel files 
install.packages("pheatmap")

# loading the libraries required 

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)
library(pheatmap)
library(ggrepel)  # for better text labels


allowWGCNAThreads() # allow multi threading ( optional)

# set the working directory 
setwd("C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\CoExpNetwork")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# 1. Fetch Data -------------------------------
data_genes <- read_excel("leaf_husk_dataset.xlsx", sheet = 1)  

data_genes
summary(data_genes)

data_long <- data_genes %>%
  select(-c(Signature, `Foliar profile`, `Husk profile`, EntrezGene)) %>%
  column_to_rownames("Accession") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Gene", values_to = "Expression")

ggplot(data_long, aes(x = Sample, y = Expression)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Gene Expression Distribution per Sample",
       x = "Sample", y = "Expression Level")

# get metadata
# no meta data?

# subseting the data to extract only the required columns 
filtered_data <- data_genes[,c(1:21)]

# convert the gene names to a uniform accession - TODO

# Convert all columns except the first to numeric
data_numeric2 <- filtered_data
data_numeric2[-1] <- lapply(filtered_data[-1], as.numeric)

data_numeric2 # columns are samples and the rows are the genes 
str(data_numeric2)


# 2. QC - Outlier detection -----------------------------------
# detect outlier genes - to exclude them from the analysis

# Transpose the data for WGCNA

# Separate the "Accession" column and transpose the numeric data only
accession_column <- data_numeric2$Accession  # Store the Accession column separately
numeric_data <- data_numeric2[-1]            # Remove the Accession column for numeric transposition

# Transpose the numeric data
data_transposed <- t(numeric_data)
data_transposed

# Assign row names based on the column names of the original data (e.g., "FP_1", "FP_2", etc.)
rownames(data_transposed) <- colnames(numeric_data)

# Assign the 'Accession' column as column names of the transposed data
colnames(data_transposed) <- accession_column

# Check the structure of data to ensure it's numeric
str(data_transposed)

# can change parameters as required, we re working with the defaults here
gsg <- goodSamplesGenes(data_transposed) # rows should be samples, columns should be genes - exact reverse of what we have - so we transpose 
summary(gsg)
gsg$allOK # if this is true, that means all the samples and vectors have no outliers 

# False - hence there are outliers

table(gsg$goodGenes) # counting the outlier genes ( can change the parameters here)
# 5585 genes are outliers

table(gsg$goodSamples) # counting all the outlier samples 
# all samples passes the outlier test 

# remove the genes that are detected as outliers 

dim(data_transposed)  # Number of rows and columns in the data
length(gsg$goodGenes) # Length of the logical vector

data_transposed <- data_transposed[, gsg$goodGenes == TRUE]
# cleaned dataset is called data_transposed from here onwards ********************

#===================================================================================
# detect outlier samples - hierarchical clustering - method 1 
htree <- hclust(dist(data_transposed), method = "average")

# Before we continue with network construction and module detection, we visualize how the clinical traits
# relate to the sample dendrogram.
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(traits, signed = FALSE, naColor = "grey")

# Save to high-resolution image (300 DPI)
png("sample_dendrogram_trait_heatmap.png", width = 3000, height = 2000, res = 300)

# # Plot the sample dendrogram and the colors underneath.
# plotDendroAndColors(htree,traitColors,
#                     groupLabels = names(traits),
#                     main = "Sample dendrogram and trait heatmap")

# Plot dendrogram and trait heatmap
plotDendroAndColors(
  dendro = hclust(dist(data_transposed), method = "average"),
  colors = traitColors,
  groupLabels = names(traits),
  main = "Sample Dendrogram and Trait Heatmap",
  cex.colorLabels = 0.8,       # Reduce trait label size
  cex.dendroLabels = 0.6,      # Reduce sample name size
  marAll = c(6, 10, 4, 2),     # Extra left margin for long trait names
  addGuide = TRUE,
  guideHang = 0.05
)

dev.off()

# =======================================================================================

# PCA - method 2 
pca <- prcomp(data_transposed) # perform a PCA 
pca.dat <- pca$x # assign the PCA data to a variable 

# finding the variance shown by each PC 
pca.var <- pca$sdev^2 # square of std dev
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2) # round to the first two decimal places 

pca.dat <- as.data.frame(pca.dat) # covert to dataframe 

# Create the PCA plot
p <- ggplot(pca.dat, aes(x = PC1, y = PC2)) +
  geom_point(color = "#0073C2FF", size = 3) +
  geom_text_repel(aes(label = rownames(pca.dat)), size = 4, max.overlaps = 10) +
  labs(
    title = "PCA Plot of Transposed Data",
    x = paste0("PC1: ", pca.var.percent[1], " %"),
    y = paste0("PC2: ", pca.var.percent[2], " %")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Display the plot
print(p)

# Save the plot in 300 dpi
ggsave("pca_plot_highres.png", plot = p, dpi = 300, width = 8, height = 6)

# outliers identified ?

# NOTE: If the samples have batch effects ( when there are many replicates), correct them before moving ahead
# TODO - do i remove the oddly clustered samples?

# exclude outlier samples - no outliers 

# 3. Normalization ---------------------------------------------------------------------

# if we have FPKM data or RPKM data, we can use them after log transforming them 
# Typically, RPKM (Reads Per Kilobase of transcript per Million mapped reads) data
# is log-transformed using the log2(x + 1) method to avoid issues with zeros in the data.

#Correct Data Structure for Normalization:
#   Rows: Genes (e.g., "GRMZM5G859187", "GRMZM5G872747", etc.)
#   Columns: Samples (e.g., "FP_1", "FP_2", etc.)

s <- t(data_transposed)

# Apply log2(x + 1) transformation to avoid log(0)
log_rpkm <- log2(t(data_transposed) + 1)

# Check the structure of the log-transformed data
str(log_rpkm)

# Optionally, transpose the log-transformed data if needed for further analysis
log_rpkm_transposed <- t(log_rpkm)

# find the zero variance genes - these could be non expressed genes or house keeping genes 
zero_variance_genes <- apply(log_rpkm_transposed, 2, var) == 0 
sum(zero_variance_genes)
# CHECK WHETHER WE HAVE TO REMOVE THESE - removed

# set the column names as the gene accessions again 
#colnames(log_rpkm_transposed) = accession_column

# this data is used for further analysis 
# to process further in the WGCNA package we require the columns to be the gene IDs and the rows to be the sample names 

# 4. Network construction -------------------------------
# Choose a set of soft - thesholding power

# function in WGCNA that helps us pick the right power that helps us build a network of scale free topology 
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
power


# Call the network topology analysis function 
sft <- pickSoftThreshold(log_rpkm_transposed,
                         powerVector = power,
                         networkType = "signed", # needs a signed network
                         verbose = 5)

# metrics calculated by this function is used to decide the power 
sft.data <- sft$fitIndices

# maximum R^2 value and the min mean.k is used to pick the best power ( scale free topology)
# pick the power 14, which is the lowest power for which the scale-free topology fit
# index curve flattens out upon reaching a high value (in this case, roughly 0.90).
# visualization to pick a power 
names (sft.data)

# Plot A: Scale-Free Topology Fit
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(color = "#0073C2FF", size = 3) +
  geom_text(nudge_y = 0.05, size = 3) +
  geom_hline(yintercept = 0.9, color = 'red', linetype = "dashed", linewidth = 1) +
  labs(
    title = "(A) Scale-Free Topology Model Fit",
    x = "Soft-Thresholding Power",
    y = "Signed R² (Scale-Free Topology Fit)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Plot B: Mean Connectivity
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point(color = "#EFC000FF", size = 3) +
  geom_text(nudge_y = 0.05, size = 3) +
  labs(
    title = "(B) Mean Connectivity by Power",
    x = "Soft-Thresholding Power",
    y = "Mean Connectivity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Combine and save the plots in high resolution
combined_plot <- grid.arrange(a1, a2, ncol = 2)

# Save as PNG with 300 dpi
ggsave("soft_thresholding_plots.png", plot = combined_plot,
       width = 12, height = 6, dpi = 300)


# we choose 14 as the soft threshold 

# create an adjacency matrix and applying soft threshold to create a weighted correlation matrix 
# using this, we calculate proximity measures - topological overlap matrix 
# used to perform hierarchical clustering 
# all this is done using one WGCNA function
# 1. pre-cluster genes into blocks ( number pre determined by the user)
# 2. full network analysis on each of these blocks separately 
# 3. the modules who are similar and have a high correlation are merged 
# 4. this method is fast and less computational intensive 

# convert matrix to numeric 
log_rpkm_transposed[] <- sapply(log_rpkm_transposed, as.numeric)

soft_power <- 12 # was 14
temp_cor <- cor # correlation function is assigned to a temp variable 
cor <- WGCNA::cor
# done to prevent errors 

#log-transformed RPKM data is stored in `log_rpkm`
num_genes <- nrow(log_rpkm)
num_genes
# 34071

sum(is.na(log_rpkm_transposed))  # Check for missing values


memory.limit(size = 64000)  # Allocate full 64GB.

# memory estimate w.r.t. block size
bwnet <- blockwiseModules(log_rpkm_transposed,
                          maxBlockSize = 35000,  # reduce block size for better memory management
                          TOMtype = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25, # threshold we want to merge the similar modules 
                          numericLabels = FALSE, # module eigen gene names will be the names of the colors instead of the numbers,
                          randomSeed = 1234, # for reproducibility,
                          verbose = 3,
                          checkMissingData = TRUE,
) 
#default is 5000
# 4G - up to 8000 - 10000
# 16G = 20 000 
# 32G = 30 000 
# in this case , we wanna process all the genes in one block that being splited into multiple blocks 

# we used default parameters here, but we need to read up and use the best parameters for the data set ( from the manual)

cor <- temp_cor # re assign the correlation function to the original correlation function 

# 5. Module Eigen genes ----------------------------------------------------------------------------------------
module_eigengenes <- bwnet$MEs # names of the module eigen genes are names of the colors.module eigengenes (bwnet$MEs) are calculated across all blocks

# Print out a preview 
head(module_eigengenes)

# get number of genes for each module - for all blocks 
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath 
# To plot the dendrograms for both blocks, you need to iterate over the blocks
# and plot each dendrogram separately with the corresponding colors for each block

# plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
#                     c("ünmerged", "merged"),
#                     dendroLabels = FALSE,
#                     addGuide = TRUE,
#                     hang = 0.03,
#                     guideHang = 0.05)
# 
# length(bwnet$dendrograms[[1]]$order)  # Number of objects in the dendrogram
# length(bwnet$unmergedColors)  # Length of unmerged colors
# length(bwnet$colors)  # Length of merged colors


# Iterate through all blocks and save dendrograms for each block as PNG files with 300 DPI
# Loop through all blocks and save high-quality dendrograms
for (block in 1:length(bwnet$dendrograms)) {
  
  # Define the filename
  filename <- paste0("Block_", block, "_GeneDendrogram_ThesisReady.png")
  
  # Extract color info for current block
  colors_to_plot <- cbind(
    bwnet$unmergedColors[bwnet$blockGenes[[block]]],
    bwnet$colors[bwnet$blockGenes[[block]]]
  )
  
  # Save plot in 300 dpi with wider margins for thesis format
  png(filename, width = 10, height = 6, units = "in", res = 300)
  
  # Plot with improved settings
  plotDendroAndColors(
    dendro = bwnet$dendrograms[[block]],
    colors = colors_to_plot,
    groupLabels = c("Before Merging", "After Merging"),
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    guideHang = 0.05,
    main = paste("Gene Clustering Dendrogram — Block", block),
    cex.main = 1.4,
    cex.rowText = 0.5
  )
  
  dev.off()
}

# the similar modules have been merged in the merged section 

# grey module - all genes that doesn't fall into other modules were assigned to the grey module 

# Relate modules to traits=================================================================
# module trait associations

# we need to turn the foliar leaf data and the husk data into binary data 

# Ensure row names are accessible
rownames(log_rpkm_transposed)

# Create the 'foliar' column: 1 if row name starts with "F", 0 otherwise
foliar <- ifelse(grepl("^F", rownames(log_rpkm_transposed)), 1, 0)

# Create the 'husk' column: 1 if row name starts with "H", 0 otherwise
husk <- ifelse(grepl("^H", rownames(log_rpkm_transposed)), 1, 0)

# Combine 'foliar' and 'husk' into a new data frame 'traits' with original row names
traits <- data.frame(foliar = foliar, husk = husk, row.names = rownames(log_rpkm_transposed))

# traits$Veins_initiated = c("Mid-vein","Mid-vein",
#                     "Mid-vein Laterals","Mid-vein Laterals",
#                     "Mid-vein Laterals Intermediates","Mid-vein Laterals Intermediates",
#                     "Mid-vein Laterals Intermediates", "Mid-vein Laterals Intermediates",
#                     "Mid-vein Laterals Intermediates","Mid-vein Laterals Intermediates",
#                     "Mid-vein","Mid-vein",
#                     "Mid-vein Laterals", "Mid-vein Laterals", 
#                     "Mid-vein Laterals", "Mid-vein Laterals", 
#                     "Mid-vein Laterals", "Mid-vein Laterals", 
#                     "Mid-vein Laterals", "Mid-vein Laterals") #- DEAL WITH THIS LATER _ TODO

traits$M_cells_between_veins = c(NA,NA,4,5,1,2, 2,2,2, 2, NA, NA, 5,5, 10, 10, 11,13,14,16) # this data was changed manually as necessary as metadata wasnt available 
traits$BS_cell_size = c(NA, NA, NA, NA, 1,1, 3,3, 3,3,NA,NA,1,1,1,1, 2,2,2,2) # small - 1, medium 2, large - 3
traits$BS_plastid_size = c(NA, NA, NA,NA,NA, NA, 1,1, 3,3, NA, NA, NA, NA, NA, NA, 1,1,1,1)
traits <- cbind(
  traits,
  data.frame(
    FP = c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    FP34 = c(0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    FP5 = c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    FI = c(0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    FE = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    HP = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    HP34 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0),
    HP5 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
    HI = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0),
    HE = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1),
    non_mature_foliar = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    mature_foliar = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    mature_husk = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
    primordia_husk = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0),
    primordia_foliar = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 , 0 , 0 , 0 ,0 , 0 , 0 , 0 , 0)
  )
)

# View the first few rows to confirm
head(traits)

# organ wise comparison - do later @TODO

# Define the number of genes and samples 
nSamples = nrow(log_rpkm_transposed)
nGenes = ncol(log_rpkm_transposed)

# Finding the module trait correlations 
module.trait.corr = cor( module_eigengenes, traits, use='p') # using Pearson correlation
module.trait.corr.pvals = corPvalueStudent(module.trait.corr, nSamples) 

# the p values helps identify which are the modules that are associated with the Foliar leaves or husk

# Create significance stars matrix
signif_stars <- ifelse(module.trait.corr.pvals < 0.001, "***",
                       ifelse(module.trait.corr.pvals < 0.01, "**",
                              ifelse(module.trait.corr.pvals < 0.05, "*", "")))

# Create combined display matrix (correlation + p value + stars)
display_matrix <- matrix(
  paste0("r=", round(module.trait.corr, 2), signif_stars, "\n", "p=", signif(module.trait.corr.pvals, 2)),
  nrow = nrow(module.trait.corr),
  dimnames = dimnames(module.trait.corr)
)

# Set up color gradient
color_gradient <- colorRampPalette(c("blue", "white", "red"))(100)

# Calculate dynamic cell sizes based on matrix dimensions
n_rows <- nrow(module.trait.corr)
n_cols <- ncol(module.trait.corr)
cell_height <- max(10, 400/n_rows)  # Minimum 10, scales with rows
cell_width <- max(30, 600/n_cols)   # Minimum 30, scales with columns

png("correlation_heatmap_new.png", 
    width = 14,         # was 12
    height = 14,        # was 8
    units = "in",       
    res = 300) 

pheatmap(module.trait.corr,
         display_numbers = display_matrix,
         number_color = "black",
         fontsize_number = 8,
         color = color_gradient,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Module-Trait Relationships",
         border_color = NA,
         cellwidth = 40,
         cellheight = 22,
         fontsize = 10,
         legend = TRUE)

dev.off()

# trait data in x axis ( that is F or H)
# eigen genes in the y axis 

# the result gives the correlations of the modules with the trait of interest
# *** means high correlation

# the red ones show which modules are related with husk 

# to filter out the genes belonging to green module( highest correlation with husk presence)
# green_module_genes <- module.gene.mapping %>%
#   filter(bwnet$colors == 'green') %>%
#   rownames()
# 
# magenta_module_genes <- module.gene.mapping %>%
#   filter(bwnet$colors == 'magenta') %>%
#   rownames()
# 
# salmon_module_genes <- module.gene.mapping %>%
#   filter(bwnet$colors == 'salmon') %>%
#   rownames()
# 
# blue_module_genes <- module.gene.mapping %>%
#   filter(bwnet$colors == 'blue') %>%
#   rownames()
# 
# yellow_module_genes <- module.gene.mapping %>%
#   filter(bwnet$colors == 'yellow') %>%
#   rownames()
# 
# we want to find the association between a trait and a module eigen gene, not the correlation
# it has to be statistically significant association ( trait vs all- whether there is a significant difference)


# 6B. Intra-modular analysis: Identifying driver genes ---------------


# Calculate the module membership and the associated p-values

# The module membership/intra-modular connectivity is calculated as the correlation of the module eigen gene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- as.data.frame(cor(log_rpkm_transposed,module_eigengenes, use = 'p'))
module.membership.measure.pvals <- as.data.frame(corPvalueStudent(as.matrix(module.membership.measure), nSamples))

# what are the genes having high module membership measures in the modules of interest 
module.membership.measure.pvals[1:10,1:10]

# # Using the module membership measures you can identify genes with high module membership in interesting modules.
# # Modules of interest
# modules_of_interest <- c("MEgreen", "MEbrown", "MEturquoise", "MEred", "MEgreenyellow", "MEsalmon", "MEblack")
# 
# # Initialize a result list to store the top genes
# top_10_genes_per_module <- list()
# 
# # Loop through each module
# for (module in modules_of_interest) {
#   
#   print(module)  # Current module
#   print(head(module_pvals))  # P-values for the current module
#   
#   if (!module %in% rownames(module.membership.measure.pvals)) {
#     stop(paste("Module", module, "not found in the p-values matrix."))
#   }
#   
#   # Get p-values for the module
#   module_pvals <- module.membership.measure.pvals[module, ]
#   
#   # Rank genes by p-value in ascending order Sorts genes by p-values in ascending order (smaller p-values are more significant).
#   ranked_genes <- sort(module_pvals, decreasing = FALSE)
#   
#   # Get the top 10 genes
#   top_10_genes <- names(ranked_genes)[1:10]
#   
#   # Store the genes and their p-values
#   top_10_genes_per_module[[module]] <- data.frame(
#     Gene = top_10_genes,
#     P_Value = ranked_genes[1:10]
#   )
# }
# 
# # View the results for each module
# top_10_genes_per_module
# 
# # Green module
# writeLines(top_10_genes_per_module$MEgreenyellow$Gene, "greenyellow_sig_genes_names.txt")
# top_10_genes_per_module$MEgreenyellow$Gene
# 
# # green
# # [1] "GRMZM2G066902" "GRMZM2G481452" "GRMZM2G092867" "GRMZM2G159013" "GRMZM2G147882"
# # [6] "GRMZM2G153877" "GRMZM2G160064" "GRMZM2G142390" "GRMZM5G846548" "GRMZM2G074472"
# 
# # > top_10_genes_per_module$MEgreenyellow$Gene
# # [1] "GRMZM5G800535" "GRMZM2G174730" "GRMZM2G410978" "GRMZM2G097468" "GRMZM2G169584" "GRMZM2G092169"
# # [7] "GRMZM2G320705" "GRMZM2G096331" "GRMZM2G063309" "GRMZM2G308597"

#====================================================================================================
# Calculate the gene significance and associated p-values ------------------------------------------

# YOU CAN CHANGE THE TRAIT HERE -  CHANGE 1 
head(traits)

gene.signf.corr <- as.data.frame(cor(log_rpkm_transposed, traits$mature_foliar , use = 'p'))

gene.signf.corr.pvals <- as.data.frame(corPvalueStudent(as.matrix(gene.signf.corr), nSamples))


#Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
# membership in interesting modules. As an example, we look at the brown module that has the highest association
# with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:

# print(unique(bwnet$colors))  # Check if "grey" exists
# sum(bwnet$colors == "green")

# module of interest
module = "turquoise" # CHANGE 2 
# extracting the color name of the modules
modNames = substring(names(module_eigengenes),3) # removing ME part
# to find which gene belongs to which module 
module.gene.mapping <- as.data.frame(bwnet$colors)
column = match(module, modNames)

# Debugging for grey module
# print(sum(moduleGenes))  # Should match the earlier 5938
# print(names(module_eigengenes))  # Confirm what the names look like
# print(modNames)  # Extracted module names
# print(module)    # The module of interest ("grey")
# print(column)    # Index of "grey" in modNames (should not be NA)
# print(unique(moduleColors))  # Check for "grey"
# WHY DID GREY DISAPPEAR?


moduleGenes = module.gene.mapping == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(module.membership.measure[moduleGenes, column]),
                   abs(gene.signf.corr[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for mature foliar stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
# CHANGE 3
# The reader is encouraged to try this code with other significance trait/module correlation

# gene Ids
names(as.data.frame(log_rpkm_transposed))
#names(as.data.frame(log_rpkm_transposed))[module.gene.mapping=="yellowgreen"]

#debugging
table(bwnet$colors)


# retrieving information 
# create a data frame to hold the information of the genes from the module of interest 
# The modules will be ordered by their significance for weight, with the most significant ones to the left.

# Create the starting data frame
geneInfo0 = data.frame(geneAccession = colnames(log_rpkm_transposed),
                       moduleColor = module.gene.mapping$`bwnet$colors`,
                       GeneSignificance = gene.signf.corr$V1,
                       GeneSignificancePvals = gene.signf.corr.pvals$V1)

# str(log_rpkm_transposed)
# str(module.gene.mapping)
# str(gene.signf.corr)
# str(gene.signf.corr.pvals)

# Order modules by their significance for weight

# CHANGE 4 HERE
modOrder = order( -abs(cor(module_eigengenes, traits$mature_husk, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(module.membership.measure))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, module.membership.measure[, modOrder[mod]],
                         module.membership.measure.pvals[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GeneSignificance));
geneInfo = geneInfo0[geneOrder, ]


# Write the data to an Excel file CHANGE 5 
write.xlsx(geneInfo, file = "geneInfo_maturehusk_darkgreen.xlsx", rowNames = FALSE)


# Exporting the network to Cytoscape 
modules = c("green", "turquoise")

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(log_rpkm_transposed, power = soft_power)
# Select module probes
probes = names(as.data.frame(log_rpkm_transposed))
moduleColors = module.gene.mapping$`bwnet$colors`
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

# ===========================================================
# set a proper threshold for this -TODO
foliar_sig_genes <- gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)


foliar_sig_genes_names <- rownames(foliar_sig_genes)
# [1] "GRMZM2G087095" "GRMZM2G010494" "GRMZM2G145554" "GRMZM2G060319" "GRMZM2G156632" "GRMZM2G132055"
# [7] "GRMZM2G700200" "GRMZM2G119465" "GRMZM2G143640" "GRMZM2G470882" "GRMZM2G030123" "GRMZM2G412436"
# [13] "GRMZM2G002976" "GRMZM2G129114" "GRMZM2G180324" "GRMZM2G026783" "GRMZM2G063262" "GRMZM2G125411"
# [19] "GRMZM2G062554" "GRMZM5G870592" "GRMZM2G018070" "GRMZM2G154093" "GRMZM2G152417" "GRMZM5G880435"
# [25] "GRMZM2G092867


# Save foliar_sig_genes_names to a text file
writeLines(foliar_sig_genes_names, "foliar_sig_genes_names.txt")

# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

# No matches found - maize GDB
# GRMZM2G002976 GRMZM2G087095 GRMZM2G700200 GRMZM2G132055 GRMZM2G030123

# NCBI
# Id=GRMZM2G132055: No items found.- plant ensemble entry found but no gene ig
# Id=GRMZM5G880435: No items found.- plant ensembl and NCBI protein entry found - manually added last to the matched genes ids file 


# STRING didnt find 
# 100192757
# 100272664
# 100273399
# 100281975
# 103642531
# ONM05113

#=======================================================================================================
# Mapping the gene names to their string IDS ========================================================

# paste the list of genes from the file to NCBI, download the result 
# manipulate the below code to extract the gene IDS 

# File paths
# insert file name here 
mapping_file <- "greenyellow_sig_genes_ncbi_output.txt"
# insert file name here
genes_of_interest_file <- "greenyellow_sig_genes_names.txt" 

# Load the data
mapping_data <- read.delim(mapping_file, header = TRUE, sep = "\t")
genes_of_interest <- readLines(genes_of_interest_file)

# Unnest Aliases for easier matching
unnested_data <- mapping_data %>%
  separate_rows(Aliases, sep = ",\\s*")

# Filter the unnested data
matched_genes <- unnested_data %>%
  filter(Aliases %in% genes_of_interest)

# Output the matched Gene IDs
matched_gene_ids <- matched_genes$GeneID

# Save Gene IDs to a plain text file, one per line - change the file name here as necessary 
writeLines(as.character(matched_gene_ids), "matched_gene_ids_greenyellow_string_input.txt")

# Print matched Gene IDs
print(matched_gene_ids)


# Module preservation analysis between c3 and c4 sample classification 
setLabels = c("c3_likeness", "c4_likeness")
c3_likeness_data = data.frame(log_rpkm_transposed[-(7:10), ]) # extract all samples except the FI and FE samples
c4_likeness_data = data.frame(log_rpkm_transposed[(7:10), ]) # extract the FI and FE samples as they were clustered together 

multiExpr=list(c3_likeness = list(data=c3_likeness_data),
               c4_likeness = list(data=c4_likeness_data))

multiColor = list(module.gene.mapping)

# The number of permutations drives the computation time
# of the module preservation function. For a publication use 200 permutations.
set.seed(1)
system.time({
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1, nPermutations = 200,
                          randomSeed = 1, quickCor = 0, verbose = 3)
})
